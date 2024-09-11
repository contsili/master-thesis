%% Simulation @ MN
% This code uses jackknifing (like the experiment). 
% Instead I can use Monte-Carlo simulation (see: dippos2_MN.m)

%% this is the code I run in the local PC

addpath('H:/common/matlab/fieldtrip')

ft_defaults



%% Same time points for OPM and SQUID (0.015*1200 +1 = 19 time points)

%% OPM

load M:\Documents\recordings\opm-MN\results\coregistration\hfc\'mat files'\headmodel_sphere.mat
load M:\Documents\recordings\opm-MN\results\dipmom_analytical\mom135.mat
load M:\Documents\recordings\opm-MN\results\dipmom_analytical\pos135.mat
load M:\Documents\recordings\opm-MN\results\jackknife\hfc\fieldlinebeta2_head_fixedsensors.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\dpu_squid.mat
load M:\Documents\recordings\opm-MN\results\jackknife\hfc\data_meg_hfc_rejectvisual_146trials.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\squid_noise_floor.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\stim_clean_146trials.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\'146 trials for OPM - no R504'\append_tl_jk1.mat

load C:\Users\kontsi\Desktop\omitted_withoutR504.mat

%% Prepare the data

% % get rid of R504 from run 3
% 
% for i=1:6
%     cfg        = [];
%     cfg.method = 'summary';
%     data_meg_hfc_clean(i) = ft_rejectvisual(cfg, data_meg_hfc_clean(i)); 
% end


%% same spatial resolution

for j=1:112  
    fieldlinebeta2_head_fixedsensors.chanpos(16+j,:)=fieldlinebeta2_head_fixedsensors.chanpos(j,:);
    fieldlinebeta2_head_fixedsensors.coilpos(16+j,:)=fieldlinebeta2_head_fixedsensors.coilpos(j,:);
    fieldlinebeta2_head_fixedsensors.coilori(16+j,:)=fieldlinebeta2_head_fixedsensors.coilori(j,:);
    fieldlinebeta2_head_fixedsensors.chanori(16+j,:)=fieldlinebeta2_head_fixedsensors.chanori(j,:);
end
%% Rename sensor that re-appear between runs

% I had 8 stable sensors between runs. I will rename them
for i = 2:6     
    index = find(ismember(data_meg_hfc_clean(i).label, {'R101_bz', 'L108_bz','R503_bz','L503_bz','L507_bz','R507_bz','R212_bz','L212_bz'}));
    for k = 1:length(index)
        data_meg_hfc_clean(i).label{index(k)} = strcat(data_meg_hfc_clean(i).label{index(k)}, '_', num2str(i));
    end    
end


% rename 'R114_bz', 'L504_bz'
for i=5:6
index = find(ismember(data_meg_hfc_clean(i).label, {'R114_bz', 'L504_bz', 'R504_bz'}));

    for k = 1:length(index)
        data_meg_hfc_clean(i).label{index(k)} = strcat(data_meg_hfc_clean(i).label{index(k)}, '_', num2str(i));
    end 
end


% rename L214_bz
for i=3
index = find(ismember(data_meg_hfc_clean(i).label, {'L214_bz'}));

    for k = 1:length(index)
        data_meg_hfc_clean(i).label{index(k)} = strcat(data_meg_hfc_clean(i).label{index(k)}, '_', num2str(i));
    end 
end



%% After I checked all the channels that I know that they are duplicate, now I will see which channels still remain duplicate
duplicate_channels = {};

for i = 1:6
    % Loop through each channel in the current run
    for j = 1:length(data_meg_hfc_clean(i).label)
        current_channel = data_meg_hfc_clean(i).label{j};
        % Compare the current channel with channels in other runs
        for k = (i + 1):6
            % Check if the channel name exists in the other runs
            if any(strcmp(current_channel, data_meg_hfc_clean(k).label))
                % If found duplicate, add it to the duplicate_channels cell array
                duplicate_channels{end + 1} = {current_channel, ['Run ', num2str(i), ' and Run ', num2str(k)]};
            end
        end
    end
end

% Display the duplicate channels
if ~isempty(duplicate_channels)
    disp('Duplicate Channels:');
    disp('Channel Name       :   Runs');
    for m = 1:length(duplicate_channels)
        disp([duplicate_channels{m}{1}, '   :   ', duplicate_channels{m}{2}]);
    end
else
    disp('No duplicate channels found.');
end


%% append data

cfg         = [];
cfg.keepsampleinfo='no';
append_data = ft_appenddata(cfg, data_meg_hfc_clean(1), data_meg_hfc_clean(2), data_meg_hfc_clean(3), data_meg_hfc_clean(4), data_meg_hfc_clean(5), data_meg_hfc_clean(6)); 

append_data.grad = fieldlinebeta2_head_fixedsensors; % QUESTION: Is is correct? Do I need to create a grad structure with 182 sensors as my sensors in append_data?

% NOTE: I have not deleted  cfg.channel     = {'all', '-L101_bz',
% '-L214_bz'}; that were not coregistered well. I think that
% fieldlinebeta2_head_fixedsensors will take care of that

%% I need to make 176 append_tl_jk.trial = L*q + add stochastic noise on my own. So basically get a 1 q and then ft_prepare_leadfield 176 times
% If I find that the relation of dpu and nr_sens is 1/sqrt() and see how it
% depends on the noise then we could use an analytical solution for our
% whole brain simulations

n = length(data_meg_hfc_clean(i).trial);

nr_sens = [16, 32, 48, 64, 80, 96, 112, 128];

% sensor noise
sensor_noise_squid = 5 * 10^-15;
% sensor_noise_squid = squid_noise_floor;
sensor_noise_opm =  3*sensor_noise_squid;

% % fit the dipmom for every time point
% cfg            = [];
% cfg.channel    = {'all', '-L101_bz', '-L214_bz'}; % delete the channels that appeared closer to the brain during the co-registration
% append_tl  = ft_timelockanalysis(cfg, append_data);
% 
% cfg            = [];
% cfg.latency    = [0.035 0.05];
% cfg.unit       = 'cm';
% cfg.gridsearch = 'yes';
% cfg.nonlinear  = 'yes'; 
% cfg.grad       = fieldlinebeta2_head_fixedsensors;                      
% cfg.headmodel  = headmodel_sphere;
% stim_dpN20_opm = ft_dipolefitting(cfg, append_tl);
% 
% mom135_146trials = mean(stim_dpN20_opm.dip.mom, 2);

% compute leadfield
cfg                 = [];      
cfg.sourcemodel.pos = pos135;
cfg.grad            = fieldlinebeta2_head_fixedsensors;
cfg.headmodel       = headmodel_sphere;
leadfield           = ft_prepare_leadfield(cfg);

% generate simulated data with sensor_noise_opm = 3*sensor_noise_squid
new_data       = append_data;
new_data.label = fieldlinebeta2_head_fixedsensors.label;
new_data.fsample = 1200;
new_data.time = stim_clean.time;
for k = 1:146    
    new_data.trial{1,k} = leadfield.leadfield{1,1} * mom135 + sensor_noise_opm * randn(length(fieldlinebeta2_head_fixedsensors.label), 0.6*1200); % note: here I add the leadfield that is defined for 35-50 ms (when the dipolar pattern is there) to the whole trial length ((-200 ms) - 400 ms). This is not correct since the dipolar pattern appears only at 35-50 ms. However I correct that later when I use cfg.latency = [0.035 0.05];
end


%% same spatial resolution
% % If I find that the relation of dpu and nr_sens is 1/sqrt() and see how it
% % depends on the noise then we could use an analytical solution for our
% % whole brain simulations
% 
% n = length(data_meg_hfc_clean(1).trial);
% 
% nr_sens = [16, 32, 48, 64, 80, 96, 112, 128];
% 
% % sensor noise
% % sensor_noise_squid = 5 * 10^-15;
% sensor_noise_squid = squid_noise_floor;
% sensor_noise_opm =  3*sensor_noise_squid;
% 
% 
% % compute leadfield
% cfg                 = [];      
% cfg.sourcemodel.pos = pos135;
% cfg.grad            = fieldlinebeta2_head_fixedsensors;
% cfg.headmodel       = headmodel_sphere;
% leadfield           = ft_prepare_leadfield(cfg);
% 
% % generate simulated data with sensor_noise_opm = 3*sensor_noise_squid
% new_data       = append_data;
% new_data.label = fieldlinebeta2_head_fixedsensors.label;
% new_data.fsample = 1200;
% new_data.time = stim_clean.time;
% 
% % Define the number of sensors and the group size
% num_sensors = length(fieldlinebeta2_head_fixedsensors.label);
% group_size = 16;
% 
% for k = 1:146    
%     new_data.trial{1,k} = repmat(leadfield.leadfield{1,1} * mom135, 1, 0.6*1200);
%     for i = 0:group_size:128
%         % Set the RNG seed for each group of 16 sensors
%         rng(1,'twister'); % Set the RNG seed based on the group index
%         % Determine the range of sensors for this group
% 
%         sensor_range = i+1:min(i+group_size, num_sensors);
% 
% 
%         % Generate noise for this group
%         noise = sensor_noise_opm * randn(length(sensor_range), 0.6*1200);
%         % Add the noise to the corresponding sensors in new_data.trial
%         new_data.trial{1,k}(sensor_range,:) = new_data.trial{1,k}(sensor_range,:) + noise;
%     end
% end



%% topo
cfg = [];
cfg.xlim = [0.035 0.050];
cfg.layout = 'fieldlinebeta2bz_helmet.mat';
cfg.marker = 'labels';
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, new_data);

%%

%  %%% 
% % % Find the indexes of elements not ending with '_bz'
% labels = append_tl_jk1;
% indexes_not_bz = find(~endsWith(labels, '_bz'));
% 
% % Generate a list of numbers from 1 to 182 excluding 'indexes_not_bz'
% all_indexes = 1:178;
% indexes_to_permute = setdiff(all_indexes, indexes_not_bz);


for m = 2:3   
    for j = 1:n
        
        % omitted = {};      
        % omitted = append_tl_jk1(indexes_to_permute(1:nr_sens(m)))';

        cfg            = [];
        cfg.trials     = setdiff(1:n, j); % leave one trial out
        cfg.channel    = {'all', '-L101_bz', '-L214_bz'}; % delete the channels that appeared closer to the brain during the co-registration
        stim_tl_jk(j)  = ft_timelockanalysis(cfg, new_data);

        cfg            = [];
        cfg.latency    = [0.035 0.05];
        cfg.unit       = 'cm';
        cfg.gridsearch = 'no';
        cfg.dip.pos    = pos135; % location where dipole was localised with 135 OPM sensors. This is the starting point for the non-linear search
        cfg.grad       = fieldlinebeta2_head_fixedsensors;
        cfg.channel    = fieldlinebeta2_head_fixedsensors.label(1:16*m);
        % cfg.channel = omitted(1:16*m);
        cfg.headmodel  = headmodel_sphere;
        stim_dpN20_opm_jk(m,j) = ft_dipolefitting(cfg, stim_tl_jk(j));

        dippos(j,:) = stim_dpN20_opm_jk(m,j).dip.pos;

        
     
    end

    bias = (n-1)^2;
    
    jackknife_std_x(m) = sqrt(bias) * std(dippos(:,1));
    jackknife_std_y(m) = sqrt(bias) * std(dippos(:,2));
    jackknife_std_z(m) = sqrt(bias) * std(dippos(:,3));
            
    dpu(m) = sqrt((jackknife_std_x(m)^2 + jackknife_std_y(m)^2 + jackknife_std_z(m)^2)/3)


end 


%% SQUID

load M:\Documents\recordings\squid-MN\results\'demeaning & ft_rejectvisual'\headmodel_sphere.mat
ctf275 = ft_read_sens('M:\Documents\recordings\squid-MN\subj001ses002_3031000.01_20231208_01.ds', 'senstype', 'meg');
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\stim_clean_146trials.mat
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\pos275.mat
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\mom275.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\squid_noise_floor.mat

%%


% % fit the dipmom for every time point
% cfg      = [];
% stim_tl  = ft_timelockanalysis(cfg, stim_clean);
% 
% cfg            = [];
% cfg.latency    = [0.035 0.05];
% cfg.unit       = 'cm';
% cfg.gridsearch = 'yes';
% cfg.nonlinear  = 'yes'; 
% cfg.grad       = ctf275;                      
% cfg.headmodel  = headmodel_sphere;
% stim_dpN20_squid = ft_dipolefitting(cfg, stim_tl);
% 
% mom275_146trials = mean(stim_dpN20_squid.dip.mom,2);


cfg                 = [];      
cfg.sourcemodel.pos = pos275;
cfg.grad            = ctf275;
cfg.headmodel       = headmodel_sphere;
leadfield_squid           = ft_prepare_leadfield(cfg); % I do not compute leadfield = ft_prepare_leadfield(cfg, append_data) since append_data already has noise inside it

sel275 = startsWith(ctf275.label, 'M');
leadfield_squid.leadfield{1,1}=leadfield_squid.leadfield{1,1}(sel275,:);

new_data_squid       = stim_clean;
new_data_squid.label = ctf275.label(sel275);
for k = 1:146    
    new_data_squid.trial{1,k} = leadfield_squid.leadfield{1,1} * mom275 + sensor_noise_squid * randn(length(ctf275.label(sel275)), 0.6*1200);
end

%% topo
cfg = [];
cfg.xlim = [0.035 0.050];
cfg.layout = 'CTF275_helmet.mat';
cfg.marker = 'labels';
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, new_data_squid);
    
%%
n = length(stim_clean.trial);

for j = 1:n

    cfg            = [];
    cfg.trials     = setdiff(1:n, j); % leave one trial out
    stim_tl_jk_squid(j)  = ft_timelockanalysis(cfg, new_data_squid);

    cfg            = [];
    cfg.latency    = [0.035 0.05];
    cfg.unit       = 'cm';
    cfg.gridsearch = 'no';  
    cfg.dip.pos    = pos275;
    cfg.grad       = ctf275;                         
    cfg.headmodel  = headmodel_sphere;
    stim_dpN20_squid_jk(j) = ft_dipolefitting(cfg, stim_tl_jk_squid(1,j));

    dippos(j,:) = stim_dpN20_squid_jk(j).dip.pos;
               

end

%%% descriptive stats of dippos. Use bias correction = (n-1)
% Biased corrected
bias = (n-1)^2;

jackknife_std_x_squid = sqrt(bias) * std(dippos(:,1));
jackknife_std_y_squid = sqrt(bias) * std(dippos(:,2));
jackknife_std_z_squid = sqrt(bias) * std(dippos(:,3));
        
dpu_squid = sqrt((jackknife_std_x_squid^2 + jackknife_std_y_squid^2 + jackknife_std_z_squid^2)/3);


%% Plot dpu/dmu - number of sensors

%% lin-lin
% Create figure with specified size
figure 

% Plot with enhanced markers
plot(nr_sens, dpu/dpu_squid, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1.5);

% fit
power_law_model = fittype('a * x^n', 'coefficients', {'a', 'n'});
fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [mean(dpu/dpu_squid), 0.5]);
fitted_model = fit(nr_sens', (dpu/dpu_squid)' , power_law_model, fit_options);

% Access the fitted coefficients
a1 = fitted_model.a;
n1 = fitted_model.n;

tmp = linspace(min(nr_sens), max(nr_sens), 1000); % More points for a smooth curve
fitted_line1=a1.*tmp.^n1;

hold on
% Plot fitted line
h_fitted_line1 =plot(tmp, fitted_line1, '--', 'Color', 'red', 'LineWidth', 2);


hold on
legend([h_fitted_line, h_fitted_line1], ...
       {sprintf('y = %.2f * x^{%.2f}', a, n), sprintf('y = %.2f * x^{%.2f}', a1, n1)}, ...
       'Location', 'best');

xticks(nr_sens);

print('7', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution

%% log-log 
%% dpu
hold on
loglog(nr_sens, dpu/dpu_squid, 's', 'MarkerFaceColor',' red', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1.5);


% Fit a line in log-log space
log_nr_sens = log10(nr_sens);
log_ratio2 = log10(dpu/dpu_squid);
coefficients = polyfit(log_nr_sens, log_ratio2, 1); % Linear fit in log-log space
a1 = coefficients(1);
b1 = coefficients(2);

% Generate fitted line data
tmp1 = linspace(min(nr_sens), max(nr_sens), 100);
fitted_line1 = 10.^(a1 * log10(tmp1) + b1);
h_fitted_line1=loglog(tmp1, fitted_line1, '--', 'Color', 'red', 'LineWidth', 1.2);

grid on

% Add horizontal line at y = 1 (which is 10^0 in log scale)
yline(1, '--', 'LineWidth', 1.5);
xticks(nr_sens);
yticks([1*(10^-1):0.1:9*(10^-1), 1*(10^0):9*(10^0), 1*(10^1):10:2*(10^1)])

% Add labels with LaTeX interpreter
xlabel('Number of OPM sensors', 'FontSize', 14); 

% Add fitted line equation as a text annotation
legend( [h_fitted_line,h_fitted_line1], sprintf('y = %.2fx + %.2f', a, b), sprintf('y = %.2fx + %.2f', -0.46, b1), 'Location', 'northeast')

% Extend x-axis and y-axis limits slightly beyond the data range
xlim([14, 135]);
ylim([0.3, 25]);

print('13', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution

