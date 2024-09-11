% New methods:

% 1. Project OPM on scalp
% 2. colocode OPM/SQUID distance from scalp and brain
% 3. FIX: *dipmom(i,:), 
% 4. Vary number of OPMs & use biaxial, triaxial ones

%%
addpath('H:/common/matlab/fieldtrip')

ft_defaults

%%

mri          = ft_read_mri('single_subj_T1_1mm.nii');
mri.coordsys = 'spm'; 

cfg                = [];
cfg.output         = {'brain','scalp'};
cfg.scalpthreshold = 0.25; 
mri_segmented      = ft_volumesegment(cfg, mri); 

cfg                = [];
cfg.tissue         = {'scalp'};
mesh_scalp         = ft_prepare_mesh(cfg, mri_segmented);

%%
figure;
ft_plot_mesh(mesh_scalp, 'facealpha', 0.5)

%% singleshell headmodel for moving opm inwards

cfg        = [];
% cfg.tissue = 'scalp';
cfg.method = 'singleshell'; 
headmodel_shell  = ft_prepare_headmodel(cfg, mesh_scalp);

%%
figure;
ft_plot_headmodel(headmodel_shell)

%% singlesphere headmodel for the forward model

cfg        = [];
% cfg.tissue = 'scalp';
cfg.method = 'singlesphere'; 
headmodel_sphere  = ft_prepare_headmodel(cfg, mesh_scalp);

%%
figure;
ft_plot_headmodel(headmodel_sphere, 'facealpha', 0.5)

%% sourcemodel

sourcemodel = ft_read_headshape('cortex_20484.surf.gii');

%% plot mri + mesh + spherical headmodel

figure;
ft_plot_headmodel(headmodel_shell, 'facealpha', 0.5, 'axes', 'true')
ft_plot_mesh(sourcemodel, 'facecolor', 'brain', 'edgecolor', 'none', 'facealpha', 0.2)
% ft_plot_sens(fieldlinebeta2)
camlight
lighting gouraud

%% sensor array opm

load fieldlinebeta2.mat; 

translation_matrix = [
    1 0 0  0;
    0 1 0  -0.025;
    0 0 1  0.01;
    0 0 0  1
];

fieldlinebeta2 = ft_transform_geometry(translation_matrix, fieldlinebeta2);

fieldlinebeta2 = ft_convert_units(fieldlinebeta2, 'mm');

hold on
ft_plot_sens(fieldlinebeta2, 'axes', 'true')

%% project opm sensors

% Iterate through each sensor
for i = 1:size(fieldlinebeta2.coilpos, 1)
    % Get the position of the current sensor
    sensor_pos = fieldlinebeta2.coilpos(i, :);
    
    % Calculate the distance from the sensor to each vertex of the scalp mesh
    distances = vecnorm(sensor_pos - mesh_scalp.pos, 2, 2);
    
    % Find the closest vertex on the scalp mesh
    [~, closest_vertex_index] = min(distances);
    closest_vertex_pos = mesh_scalp.pos(closest_vertex_index, :);
    
    % Calculate the translation vector needed to move the sensor to touch the scalp
    translation_vector = closest_vertex_pos - sensor_pos;
    
    % Apply the translation to the current sensor position
    fieldlinebeta2.coilpos(i, :) = fieldlinebeta2.coilpos(i, :) + translation_vector;
    fieldlinebeta2.chanpos(i, :) = fieldlinebeta2.coilpos(i, :);
end


%% plot opm sensors + headmodel_shell

figure
ft_plot_headmodel(headmodel_sphere, 'facealpha', 0.5, 'facecolor','skin')
ft_plot_sens(fieldlinebeta2,'edgecolor','r') % looks good
ft_plot_mesh(sourcemodel, 'facecolor', 'brain', 'edgecolor', 'none', 'facealpha', 1)
camlight
lighting gouraud

%% sensor array squid

% For the SQUID placement we follow Iivanainen's strategy: "the SQUID array covered the cortex uniformly and
% symmetrically along the left–right axis. Thereafter, the sensor positions
% were verified to be at least 2 cm (approximate thickness of the MEG
% dewar) from the scalp."

load ctf275.mat 


% rotate the sensor array
rotation_matrix = [
    0 -1 0  0;
    1 0 0  -0.04;
    0 0 1  -0.025;
    0 0 0  1
];
ctf275 = ft_transform_geometry(rotation_matrix, ctf275);

ctf275=ft_convert_units(ctf275, 'mm');

figure;
ft_plot_sens(ctf275, 'label', 'no', 'axes', 1, 'orientation', 0, 'coilshape', 'square')

%% comparison


dipmom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex');


%% opm
nr_sens = [32, 64, 96, 128];

labels = append_tl.label;
indexes_not_bz = find(~endsWith(labels, '_bz') | endsWith(labels, 'L214_bz') | endsWith(labels, 'L101_bz'));

% Generate a list of numbers from 1 to 182 excluding 'indexes_not_bz'
all_indexes = 1:182;
indexes_to_permute = setdiff(all_indexes, indexes_not_bz);

norm_lf_opm = zeros(length(sourcemodel.pos), length(nr_sens));

for m=1:length(nr_sens)
    for i = 1: length(sourcemodel.pos)
    
        omitted = {};
        omitted = append_tl.label(indexes_to_permute(1:nr_sens(m)))';
    
        cfg                 = [];
        cfg.sourcemodel.pos = sourcemodel.pos(i,:);
        cfg.grad            = fieldlinebeta2;
        cfg.channel         = omitted;
        cfg.headmodel       = headmodel_sphere;
        leadfield           = ft_prepare_leadfield(cfg);
        
        mom         = dipmom(i,:);
        norm_mom    = norm(mom);
        
        orient      = mom ./ norm_mom;
        
        lf          = leadfield.leadfield{1,1} * orient';
        norm_lf_opm(i,m)  = vecnorm(lf);

    end
end

% It took ~2 hours

%% ctf275

norm_lf_squid = zeros(length(sourcemodel.pos), 1);

for i = 1: length(sourcemodel.pos) 

    cfg                 = [];
    cfg.sourcemodel.pos = sourcemodel.pos(i,:);
    cfg.grad            = ctf275;
    cfg.headmodel       = headmodel_sphere;
    leadfield           = ft_prepare_leadfield(cfg);
    
    sel275                   = startsWith(ctf275.label, 'M');   
    leadfield.leadfield{1,1} = leadfield.leadfield{1,1}(sel275, :);

    mom         = dipmom(i,:);
    norm_mom    = norm(mom);
    
    orient      = mom ./ norm_mom;
    
    lf          = leadfield.leadfield{1,1} * orient';
    norm_lf_squid(i,1)  = vecnorm(lf);

end




%% ratio: Sensor noise of OPMs = 3* Sensor noise of SQUIDs

% ratio = zeros(length(sourcemodel.pos), length(nr_sens));

for  m=1:length(nr_sens)
    ratio(:,m) = 3* norm_lf_squid./norm_lf_opm(:,m); % dipole moment uncertainty opm/dipole moment uncertainty squid. If <1 opm is better

    log_ratio(:,m) = log10(ratio(:,m));
end

% Ensure norm of leadfield is a column vector 
% ratio=ratio(:); 

%% ratio: Sensor noise of OPMs = Sensor noise of SQUIDs

for  m=1:length(nr_sens)
    ratio(:,m) = norm_lf_squid./norm_lf_opm(:,m); % dipole moment uncertainty opm/dipole moment uncertainty squid. If <1 opm is better

    log_ratio(:,m) = log10(ratio(:,m));
end


%% plot sensors + colorcoded mesh

for m=1:length(nr_sens)
    
    figure;
    
    hs = trisurf(sourcemodel.tri, sourcemodel.pos(:, 1), sourcemodel.pos(:, 2), sourcemodel.pos(:, 3), log_ratio(:,m), 'FaceColor', 'interp','LineStyle', 'none' ); %'EdgeColor', 'none');
    camlight
    lighting gouraud
    ft_colormap('jet')
    
    caxis([min(log_ratio(:,4)), min(log_ratio(:,4))*(-1)]);
    colorbar;
    cb = colorbar;
    title(cb, '$log_{10}(\sigma_{\hat{q}, OPM}$ / $\sigma_{\hat{q}, SQUID}$)', 'Interpreter', 'latex');
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title(sprintf('%d OPM sensors', nr_sens(m)));
    grid on;
    view([90, 20])

end


%% plot middle slide to see the effect of depth

for m=1:length(nr_sens)
    figure;
    hs = trisurf(sourcemodel.tri, sourcemodel.pos(:, 1), sourcemodel.pos(:, 2), sourcemodel.pos(:, 3), log_ratio(:, m), 'FaceColor', 'interp', 'LineStyle', 'none');
    camlight;
    lighting gouraud;
    ft_colormap('jet');
    
    caxis([min(log_ratio(:, 4)), min(log_ratio(:, 4)) * (-1)]);
    colorbar;
    cb = colorbar;
    title(cb, '$\log_{10}(\sigma_{\hat{q}, OPM}$ / $\sigma_{\hat{q}, SQUID}$)', 'Interpreter', 'latex');
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title(sprintf('%d OPM sensors', nr_sens(m)));
    grid on;
    view([90, 5]);
    
    % Adjust the view to zoom in on the part where x < 0
    xlim([min(sourcemodel.pos(:, 1)), -0.3]);
end

%% OPM/SQUID DISTANCES FROM BRAIN AND SCALP

%% colorcode SQUID distance from scalp 

sel275         = startsWith(ctf275.label, 'M');
distances = pdist2(mesh_scalp.pos, ctf275.chanpos(sel275,:));
distance_squid  =  min(distances, [],2);

figure;
hs = trisurf(mesh_scalp.tri, mesh_scalp.pos(:, 1), mesh_scalp.pos(:, 2), mesh_scalp.pos(:, 3), distance_squid, 'FaceColor', 'interp','LineStyle', 'none' ); %'EdgeColor', 'none');
camlight
lighting gouraud
ft_colormap('jet')

caxis([min(distance), 30]);
colorbar
title('SQUID distance from scalp (mm)');

%% colorcode OPM distance from scalp 

distances = pdist2(mesh_scalp.pos, fieldlinebeta2.chanpos);
distance_opm  =  min(distances, [],2);

figure;
hs = trisurf(mesh_scalp.tri, mesh_scalp.pos(:, 1), mesh_scalp.pos(:, 2), mesh_scalp.pos(:, 3), distance_opm, 'FaceColor', 'interp','LineStyle', 'none' ); %'EdgeColor', 'none');
camlight
lighting gouraud
ft_colormap('jet')

colorbar
title('OPM distance from scalp (mm)');

% CONCLUSION: the opm is correctly projected on the scalp 

%% colorcode SQUID distance from brain 

sel275         = startsWith(ctf275.label, 'M');
distances = pdist2(sourcemodel.pos, ctf275.chanpos(sel275,:));
distance_squid  =  min(distances, [],2);

figure;
hs = trisurf(sourcemodel.tri, sourcemodel.pos(:, 1), sourcemodel.pos(:, 2), sourcemodel.pos(:, 3), distance_squid, 'FaceColor', 'interp','LineStyle', 'none' ); %'EdgeColor', 'none');
camlight
lighting gouraud
ft_colormap('jet')

caxis([min(distance), 50]);
colorbar
title('SQUID distance from brain (mm)');

%% colorcode OPM distance from brain 

distances = pdist2(sourcemodel.pos, fieldlinebeta2.chanpos);
distance_opm  =  min(distances, [],2);

figure;
hs = trisurf(sourcemodel.tri, sourcemodel.pos(:, 1), sourcemodel.pos(:, 2), sourcemodel.pos(:, 3), distance_opm, 'FaceColor', 'interp','LineStyle', 'none' ); %'EdgeColor', 'none');
camlight
lighting gouraud
ft_colormap('jet')

caxis([min(distance), 50]);
colorbar
title('OPM distance from brain (mm)');

%% Comparison: distance(SQUID)/distance(OPM). For deep sources ratio decreases meaning that OPM wont perform as well

sel275         = startsWith(ctf275.label, 'M');
distances = pdist2(sourcemodel.pos, ctf275.chanpos(sel275,:));
distance_squid  =  min(distances, [],2);

distances = pdist2(sourcemodel.pos, fieldlinebeta2.chanpos);
distance_opm  =  min(distances, [],2);

ratio = distance_squid./distance_opm;

figure;
hs = trisurf(sourcemodel.tri, sourcemodel.pos(:, 1), sourcemodel.pos(:, 2), sourcemodel.pos(:, 3), ratio, 'FaceColor', 'interp','LineStyle', 'none' ); %'EdgeColor', 'none');
camlight
lighting gouraud
ft_colormap('jet')

colorbar
title('Comparison: distance(SQUID)/distance(OPM)');
view([90, 20])


%% PUBLICATION FIGURES: 

%% load necessary files
load headmodel_shell.mat
sourcemodel = ft_read_headshape('cortex_20484.surf.gii');


%% Figure 3
% Plotting the SQUID sensor positions (top) and OPM sensor positions (bottom)
% Top: Frontal and lateral views of the sensor locations for one subject.

% Define the figure for the combined plots
figure('Position', [100, 100, 1200, 800]);

% Subplot for the frontal view of SQUID sensors
subplot(2, 2, 1, 'Position', [0.05, 0.55, 0.4, 0.4]);
ft_plot_headmodel(headmodel_shell, 'facecolor', 'skin', 'edgecolor', 'none', 'facealpha', 0.3, 'axes', 'false');
ft_plot_mesh(sourcemodel, 'facecolor', 'brain', 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_sens(ctf275, 'fiducial', 'false', 'coilshape', 'square');
view([0 0]); % Frontal view
camlight;
lighting gouraud;

% Subplot for the lateral view of SQUID sensors
subplot(2, 2, 2, 'Position', [0.55, 0.55, 0.4, 0.4]);
ft_plot_headmodel(headmodel_shell, 'facecolor', 'skin', 'edgecolor', 'none', 'facealpha', 0.3, 'axes', 'false');
ft_plot_mesh(sourcemodel, 'facecolor', 'brain', 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_sens(ctf275, 'fiducial', 'false', 'coilshape', 'square');
view([90 0]); % Lateral view
camlight;
lighting gouraud;

% Subplot for the frontal view of OPM sensors
subplot(2, 2, 3, 'Position', [0.05, 0.05, 0.4, 0.4]);
ft_plot_headmodel(headmodel_shell, 'facecolor', 'skin', 'edgecolor', 'none', 'facealpha', 0.3, 'axes', 'false');
ft_plot_mesh(sourcemodel, 'facecolor', 'brain', 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_sens(fieldlinebeta2, 'fiducial', 'false', 'coilshape', 'point', 'coilsize', 20);
view([0 0]); % Frontal view
camlight;
lighting gouraud;

% Subplot for the lateral view of OPM sensors
subplot(2, 2, 4, 'Position', [0.55, 0.05, 0.4, 0.4]);
ft_plot_headmodel(headmodel_shell, 'facecolor', 'skin', 'edgecolor', 'none', 'facealpha', 0.3, 'axes', 'false');
ft_plot_mesh(sourcemodel, 'facecolor', 'brain', 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_sens(fieldlinebeta2, 'fiducial', 'false', 'coilshape', 'point', 'coilsize', 20);
view([90 0]); % Lateral view
camlight;
lighting gouraud;

print('3', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution

%% Figure 4

% Create a figure with 2 rows and 4 columns
figure;

for m = 1:length(nr_sens)
    % Plot lateral view (first row)
    subplot(2, 4, m);
    trisurf(sourcemodel.tri, sourcemodel.pos(:, 1), sourcemodel.pos(:, 2), sourcemodel.pos(:, 3), ...
            log_ratio(:, m), 'FaceColor', 'interp', 'LineStyle', 'none');
    camlight;
    lighting gouraud;
    colormap('jet');
    % caxis([min(log_ratio(:, end)), -min(log_ratio(:, end))]);
    caxis([-0.155, 0.155]);
    view([90, 20]);
    axis off;
    daspect([1 1 1]); % Set data aspect ratio
    title(sprintf('%d OPM sensors', nr_sens(m)));
    
    % Plot mid-sagittal view (second row)
    subplot(2, 4, m + length(nr_sens));
    trisurf(sourcemodel.tri, sourcemodel.pos(:, 1), sourcemodel.pos(:, 2), sourcemodel.pos(:, 3), ...
            log_ratio(:, m), 'FaceColor', 'interp', 'LineStyle', 'none');
    camlight;
    lighting gouraud;
    colormap('jet');
    % caxis([min(log_ratio(:, end)), -min(log_ratio(:, end))]);
    caxis([-0.155, 0.155]);
    view([90, 5]);
    axis off;
    daspect([1 1 1]); % Set data aspect ratio

    %  Plot mid-sagittal slice
    xlim([min(sourcemodel.pos(:, 1)), -0.3]);
end

% Add colorbar in a separate axis
h = colorbar('Position', [0.92 0.1 0.02 0.8]); % Adjust the position to fit your layout

% ticks = linspace(min(log_ratio(:, end)), -min(log_ratio(:, end)), 5); % 5 ticks
ticks = linspace(-0.155, 0.155, 5); % 5 ticks

tickLabels = arrayfun(@(x) sprintf('%.2f', 10^x), ticks, 'UniformOutput', false);
set(h, 'Ticks', ticks, 'TickLabels', tickLabels);

set(gcf, 'Position', [100, 100, 1200, 600]); % Adjust figure size to fit the aspect ratio

print('15', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution
