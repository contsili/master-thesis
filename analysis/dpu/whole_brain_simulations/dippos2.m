%% New methods: 
% 1a. Use the fieldlinebeta2 and ctf275 array as we did in dipmom5
% 1b. do not use code from NI2 course
% 2. Pick sensor noise OPM = 3 * sensor noise SQUID
% 3. MN simulation & Whole brain simulation 

%% start fieldtrip

addpath('H:/common/matlab/fieldtrip')

ft_defaults

%% load
load H:/megmethods/kontsi/Documents/recordings/simulations/dipmom/'mat files'/mri_segmented.mat
load H:/megmethods/kontsi/Documents/recordings/simulations/dipmom/'mat files'/append_tl.mat

addpath H:/megmethods/kontsi/Documents/recordings/simulations/dipmom
%%
cfg = [];
cfg.tissue = {'scalp'};
mesh_scalp = ft_prepare_mesh(cfg, mri_segmented);
%% singlesphere headmodel for the forward model
cfg = [];
cfg.method = 'singlesphere';
headmodel_sphere = ft_prepare_headmodel(cfg, mesh_scalp);

cfg = [];
cfg.method = 'singleshell';
headmodel_shell = ft_prepare_headmodel(cfg, mesh_scalp);


%% sourcemodel
% sourcemodel = ft_read_headshape('cortex_20484.surf.gii');
sourcemodel = ft_read_headshape('cortex_5124.surf.gii');
%% sensor array opm
load fieldlinebeta2.mat;
translation_matrix = [
1 0 0 0;
0 1 0 -0.025;
0 0 1 0.01;
0 0 0 1
];
fieldlinebeta2 = ft_transform_geometry(translation_matrix, fieldlinebeta2);
fieldlinebeta2 = ft_convert_units(fieldlinebeta2, 'mm');

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
%% sensor array squid
load ctf275.mat
% rotate the sensor array
rotation_matrix = [
0 -1 0 0;
1 0 0 -0.04;
0 0 1 -0.025;
0 0 0 1
];
ctf275 = ft_transform_geometry(rotation_matrix, ctf275);
ctf275=ft_convert_units(ctf275, 'mm');

%% plot
figure
ft_plot_sens(fieldlinebeta2)
hold on
ft_plot_mesh(sourcemodel, 'facealpha', 0.9)
scatter3(sourcemodel.pos(16000,1), sourcemodel.pos(16000,2), sourcemodel.pos(16000,3), 'filled', 'SizeData', 100)
quiver3(sourcemodel.pos(16000,1), sourcemodel.pos(16000,2), sourcemodel.pos(16000,3), ...
        dipmom(16000,1), dipmom(16000,2), dipmom(16000,3), 'r', 'LineWidth', 2);
axis on
xlabel('x')
ylabel('y')
camlight
lighting gouraud

%% Sensor noise
sensor_noise_squid = 5 * 10^-15;
sensor_noise_opm =  3*sensor_noise_squid;

%% comparison

dipmom = surface_normals(sourcemodel.pos, sourcemodel.tri, 'vertex');

%% opm

nr_sens = [32, 64, 96, 128];

labels = append_tl.label;
indexes_not_bz = find(~endsWith(labels, '_bz') | endsWith(labels, 'L214_bz') | endsWith(labels, 'L101_bz'));

% Generate a list of numbers from 1 to 182 excluding 'indexes_not_bz'
all_indexes = 1:182;
indexes_to_permute = setdiff(all_indexes, indexes_not_bz);


% sigma_location_opm = zeros(length(sourcemodel.pos), length(nr_sens)); % initialise vector

for m = 1:length(nr_sens)
    for i = 1
    
    omitted = {};
    omitted = append_tl.label(indexes_to_permute(1:nr_sens(m)))';

    %%%%%%%% Leadfield %%%%%%%%
    cfg = [];
    cfg.sourcemodel.pos = sourcemodel.pos(i,:);
    cfg.grad = fieldlinebeta2;
    cfg.channel = omitted;
    cfg.headmodel = headmodel_sphere;
    leadfield = ft_prepare_leadfield(cfg);
    
    mom = dipmom(i,:);
    norm_mom = norm(mom);
    orient = mom ./ norm_mom;
    lf = leadfield.leadfield{1,1} * orient';
   
    %%%%% Each sensor has each own sensor noise that does not change for
    %%%%% different m or i %%%%
    rng(1, 'twister') % produce the same set of random numbers in each iteration

    %%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%
    data = [];
    data.time=linspace(0.01,1,100);
    data.avg = repmat(lf,1,length(data.time))+ sensor_noise_opm * randn(nr_sens(m), length(data.time));   %%%%%%%% Add Gaussian noise to leadfield %%%%%%%
    data.label = leadfield.label; % It is wrong to use data.label = omitted. ft_prepare_leadfield() reshuffles the labels and it has leadfield.leadfield{1,1}(2,:) corresponds to leadfield.label(2)=R102_bz and NOT to omitted(2) = R210_bz
    data.grad = fieldlinebeta2; % Note: use data.elec for eeg and data.grad for meg
    data.dimord = 'chan_time';
    
    cfg = [];
    cfg.gridsearch = 'no';
    cfg.dip.pos = sourcemodel.pos(i,:);
    cfg.model = 'moving';
    cfg.latency = 'all';
    cfg.headmodel = headmodel_sphere;
    cfg.nonlinear = 'yes';
    cfg.numdipoles = 1;
    cfg.showcallinfo = 'no';
    
    dip = ft_dipolefitting(cfg, data); 
    


    % % way1: uses MAP_toolbox which has limited lisences in DCCN
    % from_structure_to_vector = extractfield(dip.dip, 'pos');
    % dip_location_opm = reshape(from_structure_to_vector, [3,length(data.time)])';

    % way2
    dip_location_opm = cat(1,dip.dip.pos);

    sigma_location_opm(i,m) = sqrt((std(dip_location_opm(:,1))^2 + std(dip_location_opm(:,2))^2 + std(dip_location_opm(:,3))^2)  / 3 ); % this formula is taken from Vrba (2000). "Multichannel SQUID biomagnetic systems"

    end 
end

figure;
plot(nr_sens, sigma_location_opm*10^-1)
xlabel('N')
ylabel('Standard Deviation of Dipole Position (cm)')

%% topo of 128 channels

cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.channel         = {'all', '-R105_bz'	'-R603_bz'	'-R604_bz'	'-R605_bz'	'-L101_bz'	'-L105_bz'	'-L109_bz'	'-L113_bz'	'-L206_bz'	'-L213_bz'	'-L214_bz'	'-L307_bz'	'-L308_bz'	'-L312_bz'	'-L509_bz'	'-L604_bz' }; 

data        = [];
data.avg    = dip.Vdata(:,1); % this is the same as data.avg = data.avg(:,1)
data.label  = dip.label;
data.time   = 1;
data.dimord = 'chan_time';

ft_topoplotER(cfg, data)
 
%% squid

sigma_location_squid = zeros(length(sourcemodel.pos),1); % initialise vector

for i = 1         
    %%%%%%%% Leadfield %%%%%%%%
    cfg                 = [];
    cfg.sourcemodel.pos = sourcemodel.pos(i,:);
    cfg.grad            = ctf275;
    cfg.headmodel = headmodel_sphere;
    leadfield = ft_prepare_leadfield(cfg);

    sel275 = startsWith(ctf275.label, 'M');
    leadfield.leadfield{1,1}=leadfield.leadfield{1,1}(sel275,:);
    
    mom = dipmom(i,:);
    norm_mom = norm(mom);
    orient = mom ./ norm_mom;
    lf275 = leadfield.leadfield{1,1} * orient';


    
    %%%%% Each sensor has each own sensor noise that does not change for
    %%%%% different i %%%%
    rng(1, 'twister') % produce the same set of random numbers in each iteration

    %%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%
    data = [];
    data.time=linspace(0.01,1,100);
    data.avg = repmat(lf275,1,length(data.time))+ sensor_noise_squid * randn(275, length(data.time)); 
    data.label = ctf275.label(sel275,:);
    data.grad = ctf275; % Note: use data.elec for eeg and data.grad for meg
    data.dimord = 'chan_time';
    
    cfg = [];
    cfg.gridsearch = 'no';
    cfg.dip.pos = sourcemodel.pos(i,:);
    cfg.model = 'moving';
    cfg.latency = 'all';
    cfg.headmodel = headmodel_sphere;
    cfg.nonlinear = 'yes';
    cfg.numdipoles = 1;
    cfg.showcallinfo = 'no';
    
    dip = ft_dipolefitting(cfg, data); 
    
    
    dip_location_squid = cat(1,dip.dip.pos);

    sigma_location_squid(i,:) = sqrt((std(dip_location_squid(:,1))^2 + std(dip_location_squid(:,2))^2 + std(dip_location_squid(:,3))^2)  / 3 ); % this formula is taken from Vrba (2000). "Multichannel SQUID biomagnetic systems"

end 

    
%% topo squid

figure

cfg = [];
cfg.layout = 'CTF275_helmet';

data        = [];
data.avg    = dip.Vdata(:,1);
data.label  = dip.label;
data.time   = 1;
data.dimord = 'chan_time';

ft_topoplotER(cfg, data)

%% take care of the outputs from the HFC
sigma_location_opm_final = zeros(length(sourcemodel.pos),length(nr_sens));

load sigma_location_opm1.mat
sigma_location_opm1 = sigma_location_opm;

load sigma_location_opm2.mat
sigma_location_opm1 = [sigma_location_opm1, sigma_location_opm];

for i = 1:length(sigma_location_opm1)
    sigma_location_opm_final = sigma_location_opm1{1,i}+sigma_location_opm_final;
end

sigma_location_squid_final = zeros(length(sourcemodel.pos),1);

load sigma_location_squid1.mat
sigma_location_squid1 = sigma_location_squid;

load sigma_location_squid2.mat
sigma_location_squid1 = [sigma_location_squid1, sigma_location_squid];

for i = 1:length(sigma_location_squid1)
    sigma_location_squid_final = sigma_location_squid1{1,i}+sigma_location_squid_final;
end

%% ratio

for m=1:length(nr_sens)
    ratio(:,m) = sigma_location_opm_final(:,m)./sigma_location_squid_final; % dipole position uncertainty opm/dipole position uncertainty squid. If <1 opm is better
    log_ratio(:,m) = log10(ratio(:,m));
end

%% 1D plot

% plot(nr_sens, ratio)

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
    % hold on;
    % ft_plot_sens(fieldlinebeta2, 'label', 'no', 'axes', 0, 'orientation', 0)
    %
    % hold on
    % ft_plot_sens(ctf275, 'label', 'no', 'axes', 0, 'orientation', 0, 'coilshape', 'square')
    %
    % hold on
    % ft_plot_mesh(mesh_scalp, 'facealpha', 0.3, 'facecolor','skin', 'edgecolor', 'skin')
    % camlight
    % lighting gouraud
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


%% put one sensor on top of each other
%% opm

nr_sens = [32, 64, 96, 128];

labels = append_tl.label;
indexes_not_bz = find(~endsWith(labels, '_bz') | endsWith(labels, 'L214_bz') | endsWith(labels, 'L101_bz'));

% Generate a list of numbers from 1 to 182 excluding 'indexes_not_bz'
all_indexes = 1:182;
indexes_to_permute = setdiff(all_indexes, indexes_not_bz);


% sigma_location_opm_spatialres = zeros(length(sourcemodel.pos), length(nr_sens)); % initialise vector

for m = 1:length(nr_sens)
    for i = 1:1
    
    % for simulation (one sensor on top of another)
    omitted = {};
    omitted = append_tl.label(indexes_to_permute(1:32))'; %%%% 32 sensors on top of each other %%%%%

    %%%%%%%% Leadfield %%%%%%%%
    cfg = [];
    cfg.sourcemodel.pos = sourcemodel.pos(i,:);
    cfg.grad = fieldlinebeta2;
    cfg.channel = omitted;
    cfg.headmodel = headmodel_sphere;
    leadfield = ft_prepare_leadfield(cfg);
    
    leadfield.leadfield{1,1} = m*leadfield.leadfield{1,1}; %%%%% leadfield(2 sensors on top of each other) = leadfield(1 sensor) + leadfield(1 sensor that is on top of the previous sensor) %%%%%%

    mom = dipmom(i,:);
    norm_mom = norm(mom);
    orient = mom ./ norm_mom;
    lf = leadfield.leadfield{1,1} * orient';
   
  
    %%%%% Each sensor has each own sensor noise that does not change when
    %%%%% we add another sensor on top of it %%%%%%%
    rng(1, 'twister') % produce the same set of random numbers in each iteration


    %%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%
    data = [];
    data.time=linspace(0.01,1,100);
    data.avg = repmat(lf,1,length(data.time))+ sensor_noise_opm * randn(length(omitted), length(data.time));   %%%%%%%% Add Gaussian noise to leadfield %%%%%%%
    data.label = leadfield.label;
    data.grad = fieldlinebeta2; 
    data.dimord = 'chan_time';
    
    cfg = [];
    cfg.gridsearch = 'no';
    cfg.dip.pos = sourcemodel.pos(i,:);
    cfg.model = 'moving';
    cfg.latency = 'all';
    cfg.headmodel = headmodel_sphere;
    cfg.nonlinear = 'yes';
    cfg.numdipoles = 1;
    cfg.showcallinfo = 'no';
    
    dip = ft_dipolefitting(cfg, data); 
    
    %% topos

    % cfg = [];
    % cfg.layout = 'fieldlinebeta2bz_helmet';
    % 
    % data1        = [];
    % data1.avg    = dip.Vdata(:,1); % this is the same as data.avg = data.avg(:,1)
    % data1.label  = dip.label;
    % data1.time   = 1;
    % data1.dimord = 'chan_time';
    % 
    % ft_topoplotER(cfg, data1)


    % cfg = [];
    % cfg.layout = 'fieldlinebeta2bz_helmet';
    % 
    % data1        = [];
    % data1.avg    = abs(dip.Vmodel(:,1)- dip.Vdata(:,1)); % this is the same as data.avg = data.avg(:,1)
    % data1.label  = dip.label;
    % data1.time   = 1;
    % data1.dimord = 'chan_time';
    % 
    % ft_topoplotER(cfg, data1)
    % title('residue')

    %%
    from_structure_to_vector = extractfield(dip.dip, 'pos');
    dip_location_opm = reshape(from_structure_to_vector, [3,length(data.time)])'; 

    sigma_location_opm_spatialres(i,m) = sqrt((std(dip_location_opm(:,1))^2 + std(dip_location_opm(:,2))^2 + std(dip_location_opm(:,3))^2)  / 3 ); % this formula is taken from Vrba (2000). "Multichannel SQUID biomagnetic systems"

    end 
end



%% lin-lin

figure;
plot(nr_sens, sigma_location_opm_spatialres/sigma_location_squid(1,:), 'o', 'MarkerFaceColor', 'RED');
hold on
yline(1, '--')

xlabel('N')
ylabel('$\sigma_{\hat{r}, OPM}$ / $\sigma_{\hat{r}, SQUID}$ ', 'Interpreter', 'latex', 'FontSize', 14);
title('Simulation for Dipole at random location on cortical sheet')

%% log-log 

figure;
plot(log10(nr_sens), log10(sigma_location_opm_spatialres/sigma_location_squid(1,:)), 'o','MarkerFaceColor', 'RED');
hold on
yline(0, '--')
xlabel('$\log_{10}(N)$','Interpreter', 'latex', 'FontSize', 14);
ylabel('$log_{10}(\sigma_{\hat{r}, OPM}$ / $\sigma_{\hat{r}, SQUID}$)', 'Interpreter', 'latex', 'FontSize', 14);
title('Log-log plot: Simulation for Dipole at random location on cortical sheet')

%%% fit a line to the log-log plot. The slope is -0.5
coefficients = polyfit(log10(nr_sens), log10(sigma_location_opm_spatialres/sigma_location_squid(1,:)), 1); % Linear fit (degree 1)
a = coefficients(1);
b = coefficients(2);

tmp = linspace(30,510,10);
fitted_line=a*log10(tmp)+b;
plot(log10(tmp), fitted_line, '--');

fitted_line_eq = sprintf('Fitted line: y = %.2fx + %.2f', a, b);
legend(fitted_line_eq, 'Location', 'northwest');

%% PUBLICATION FIGURES
%% Figure 5

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
    caxis([min(log_ratio(:, end)), -min(log_ratio(:, end))]);
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
    caxis([min(log_ratio(:, end)), -min(log_ratio(:, end))]);
    view([90, 5]);
    axis off;
    daspect([1 1 1]); % Set data aspect ratio

    %  Plot mid-sagittal slice
    xlim([min(sourcemodel.pos(:, 1)), -0.3]);
end

% Add colorbar in a separate axis
h = colorbar('Position', [0.92 0.1 0.02 0.8]); % Adjust the position to fit your layout

% Set colorbar ticks and labels
% min_val = min(log_ratio(:, end));
% max_val = -min(log_ratio(:, end));
% log_ticks = [log10(0.5), log10(0.6), log10(0.7), log10(0.8), log10(0.9), 0, log10(2)];
% ticks = [min_val, log_ticks, max_val]; % Combine with min_val and max_val
% tickLabels = arrayfun(@(x) sprintf('%.2f', 10^x), ticks, 'UniformOutput', false);
% set(h, 'Ticks', ticks, 'TickLabels', tickLabels);

ticks = linspace(min(log_ratio(:, end)), -min(log_ratio(:, end)), 5); % 5 ticks
tickLabels = arrayfun(@(x) sprintf('%.2f', 10^x), ticks, 'UniformOutput', false);
set(h, 'Ticks', ticks, 'TickLabels', tickLabels);

title(h, '\textbf{$\log_{10}(DPU_{OPM}$ / $DPU_{SQUID}$})', 'Interpreter', 'latex','FontSize', 8);

set(gcf, 'Position', [100, 100, 1200, 600]); % Adjust figure size to fit the aspect ratio

print('5', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution

