%% Analytical method: uses forward model with uncorrelated sensor noise

% Add FieldTrip to the path
addpath('H:\common\matlab\fieldtrip')

ft_defaults


%% SQUID

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\pos275.mat
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\mom275.mat

load M:\Documents\recordings\squid-MN\results\'demeaning & ft_rejectvisual'\headmodel_sphere.mat
ctf275 = ft_read_sens('M:\Documents\recordings\squid-MN\subj001ses002_3031000.01_20231208_01.ds', 'senstype', 'meg');

%%% norm_lf
cfg =[];
cfg.sourcemodel.pos = pos275;
cfg.grad            = ctf275; 
cfg.headmodel       = headmodel_sphere;
leadfield_squid     = ft_prepare_leadfield(cfg);

sel275 = startsWith(ctf275.label, 'M');
leadfield_squid.leadfield{1,1}=leadfield_squid.leadfield{1,1}(sel275,:);

% QUESTION: leadfield.leadfield{1,1} is a (nr of sensors)x3 matrix. How am I
% going to calculate 1 number for the norm and integrate the x,y and z orientations of the leadfield?

mom = mom275;
norm_mom = norm(mom);
orient = mom ./ norm_mom;


lf = leadfield_squid.leadfield{1,1} * orient;
norm_lf_squid = vecnorm(lf);

%% load

load M:\Documents\recordings\opm-MN\results\preprocessing\hfc\append_tl.mat

load M:\Documents\recordings\opm-MN\results\coregistration\hfc\'mat files'\headmodel_sphere.mat
load M:\Documents\recordings\opm-MN\results\dipmom_analytical\pos135.mat
load M:\Documents\recordings\opm-MN\results\dipmom_analytical\mom135.mat

load M:\Documents\recordings\opm-MN\results\jackknife\hfc\fieldlinebeta2_head_fixedsensors.mat


% load M:\Documents\recordings\opm-MN\results\dipmom_analytical\dmu_squid.mat
% load M:\Documents\recordings\opm-MN\results\dipmom_analytical\norm_lf_squid.mat

%%

% In the plot I got from the experiment 3 things vary: 
% 1. number of OPM sensors, 
% 2. sensor noise profiles of each sensor (e.g., noise_sens1 > noise_sens2)
% 3. the spatial placement of each sensor (e.g., sensor close to S1 have higher L)
% 
% TODO To get a perfect 1/sqrt(nr_sens) relation I need to control 2. and 3.

nr_sens = [16, 32, 48, 64, 80, 96, 112, 128];

%%

labels = append_tl.label;
indexes_not_bz = find(~endsWith(labels, '_bz') | endsWith(labels, 'L214_bz') | endsWith(labels, 'L101_bz'));

% Generate a list of numbers from 1 to 182 excluding 'indexes_not_bz'
all_indexes = 1:182;
indexes_to_permute = setdiff(all_indexes, indexes_not_bz);

%% same spatial resolution

for j=1:112  
    fieldlinebeta2_head_fixedsensors.chanpos(16+j,:)=fieldlinebeta2_head_fixedsensors.chanpos(j,:);
    fieldlinebeta2_head_fixedsensors.coilpos(16+j,:)=fieldlinebeta2_head_fixedsensors.coilpos(j,:);
    fieldlinebeta2_head_fixedsensors.coilori(16+j,:)=fieldlinebeta2_head_fixedsensors.coilori(j,:);
    fieldlinebeta2_head_fixedsensors.chanori(16+j,:)=fieldlinebeta2_head_fixedsensors.chanori(j,:);
end

%% CONTROL 2.: All the sensors have the same sensor noise profiles: noise_opm=3*noise_squid 


for m = 1:length(nr_sens)    
    
    omitted = {};
    omitted = append_tl.label(indexes_to_permute(1:nr_sens(m)))';

    
    %%% norm_lf
    cfg                 = [];
    cfg.sourcemodel.pos = pos135;
    cfg.grad            = fieldlinebeta2_head_fixedsensors;
    cfg.channel         = omitted;
    cfg.headmodel       = headmodel_sphere;
    leadfield           = ft_prepare_leadfield(cfg);
    
    mom         = mom135;
    norm_mom    = norm(mom);    
    orient      = mom ./ norm_mom;
    
    lf          = leadfield.leadfield{1,1} * orient;
    norm_lf(m)  = vecnorm(lf);
  
    ratio2(m)     = 3*norm_lf_squid / norm_lf(m) ;
    log_ratio2(m) = log10(ratio2(m));
    
end



%% PUBLICATION FIGURES
%% lin-lin

% Create figure with specified size
figure('Units', 'Inches', 'Position', [1, 1, 6, 4]); 

hold on
% Plot with enhanced markers
plot(nr_sens, ratio2, 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1.5);

% fit
power_law_model = fittype('a * x^n', 'coefficients', {'a', 'n'});
fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [mean(ratio2), 0.5]);
fitted_model = fit(nr_sens', ratio2' , power_law_model, fit_options);

% Access the fitted coefficients
a = fitted_model.a;
n = fitted_model.n;

tmp = linspace(min(nr_sens), max(nr_sens), 1000); % More points for a smooth curve
fitted_line=a.*tmp.^n;

hold on
% Plot fitted line
h_fitted_line =plot(tmp, fitted_line, '--', 'Color', 'blue', 'LineWidth', 2);


% Add grid and box
grid on;
box on;

% Set axis labels
xlabel('OPM sensors', 'FontSize', 14);
ylabel('\textbf{$DMU_{OPM}$ / $DMU_{SQUID}$}', 'Interpreter', 'latex','FontSize', 8);
ylim([0.5, max(ratio2)+0.2]);

hold on
yline(1, '--','LineWidth', 2)
legend(h_fitted_line , sprintf('y = %.2f * x^{%.2f}', a, n), 'Location', 'best');

print('6', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution



%% log-log

figure;
loglog(nr_sens, ratio2, 's', 'MarkerFaceColor',' blue', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on

% Fit a line in log-log space
log_nr_sens = log10(nr_sens);
log_ratio2 = log10(ratio2);
coefficients = polyfit(log_nr_sens, log_ratio2, 1); % Linear fit in log-log space
a = coefficients(1);
b = coefficients(2);

% Generate fitted line data
tmp = linspace(min(nr_sens), max(nr_sens), 100);
fitted_line = 10.^(a * log10(tmp) + b);
h_fitted_line=loglog(tmp, fitted_line, '--', 'Color', 'blue', 'LineWidth', 1.2);


