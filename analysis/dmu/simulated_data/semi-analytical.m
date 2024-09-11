%% Semi-analytical method: uses forward model with uncorrelated sensor noise, correlated brain and environmental noise

%% this is the code I run in the local PC

addpath('H:/common/matlab/fieldtrip')

ft_defaults

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

%% same spatial resolution

for j=1:96  
    fieldlinebeta2_head_fixedsensors.chanpos(32+j,:)=fieldlinebeta2_head_fixedsensors.chanpos(j,:);
    fieldlinebeta2_head_fixedsensors.coilpos(32+j,:)=fieldlinebeta2_head_fixedsensors.coilpos(j,:);
    fieldlinebeta2_head_fixedsensors.coilori(32+j,:)=fieldlinebeta2_head_fixedsensors.coilori(j,:);
    fieldlinebeta2_head_fixedsensors.chanori(32+j,:)=fieldlinebeta2_head_fixedsensors.chanori(j,:);
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

%% Same time points for OPM and SQUID (0.015*1200 +1 = 19 time points)


n = length(data_meg_hfc_clean(i).trial);

nr_sens = [32,64,96,128];

% sensor noise
sensor_noise_squid = squid_noise_floor;
sensor_noise_opm   =  3*sensor_noise_squid;

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





%% 

n = length(data_meg_hfc_clean(i).trial);

nr_sens = [32,64,96,128]; 

 %%% 
% % Find the indexes of elements not ending with '_bz'
labels = append_tl_jk1;
indexes_not_bz = find(~endsWith(labels, '_bz'));

% Generate a list of numbers from 1 to 182 excluding 'indexes_not_bz'
all_indexes = 1:174;
indexes_to_permute = setdiff(all_indexes, indexes_not_bz);
%%%  

% dmu3
mom = mom135; % data.avg (simulated signal) needs lf. lf needs orient (dipole orientation). orient needs mom (dipole moment). In an ideal simulation I generate a dipole and I know its mom. Here however, the best estimation I have is mom135 (estimation of moment from 135 OPM sensors)
norm_mom = norm(mom);
orient = mom ./ norm_mom;



for m = 1:length(nr_sens)    
  
    
    omitted = {};      
    omitted = append_tl_jk1(indexes_to_permute(1:nr_sens(m)))';

    cfg            = [];
    cfg.channel    = {'all', '-L101_bz', '-L214_bz'}; % delete the channels that appeared closer to the brain during the co-registration
    stim_tl_jk  = ft_timelockanalysis(cfg, new_data);

    cfg            = [];
    cfg.latency    = [0.035 0.05];
    cfg.unit       = 'cm';
    cfg.gridsearch = 'no';
    cfg.dip.pos    = pos135; % location where dipole was localised with 135 OPM sensors. This is the starting point for the non-linear search
    cfg.grad       = fieldlinebeta2_head_fixedsensors;
    % cfg.channel    = fieldlinebeta2_head_fixedsensors.label(1:32*m);
    cfg.channel=omitted;
    cfg.headmodel  = headmodel_sphere;
    stim_dpN20_opm_jk(m,:) = ft_dipolefitting(cfg, stim_tl_jk);

    dippos = stim_dpN20_opm_jk(m,:).dip.pos;

    % dmu3
    cfg                 = [];
    cfg.sourcemodel.pos = pos135;
    cfg.grad            = fieldlinebeta2_head_fixedsensors;
    cfg.channel         = omitted;
    cfg.headmodel       = headmodel_sphere;
    leadfield(m)        = ft_prepare_leadfield(cfg);

    dipmom(m,:) = pinv(leadfield(m).leadfield{1,1}*orient)*stim_dpN20_opm_jk(m,:).Vdata ;

     
    dmu(m) = sqrt(sum((dipmom(m,:)-norm_mom).^2)/19);

end



%% SQUID

load M:\Documents\recordings\squid-MN\results\'demeaning & ft_rejectvisual'\headmodel_sphere.mat
ctf275 = ft_read_sens('M:\Documents\recordings\squid-MN\subj001ses002_3031000.01_20231208_01.ds', 'senstype', 'meg');
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\stim_clean_146trials.mat
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\pos275.mat
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\mom275.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\squid_noise_floor.mat

%%

sensor_noise_squid = squid_noise_floor;


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

%%
n = length(stim_clean.trial);

% dmu3
mom = mom275; 
norm_mom_squid = norm(mom);
orient = mom ./ norm_mom_squid;


cfg            = [];
stim_tl_jk_squid  = ft_timelockanalysis(cfg, new_data_squid);

cfg            = [];
cfg.latency    = [0.035 0.05];
cfg.unit       = 'cm';
cfg.gridsearch = 'no';  
cfg.dip.pos    = pos275;
cfg.grad       = ctf275;                         
cfg.headmodel  = headmodel_sphere;
stim_dpN20_squid_jk = ft_dipolefitting(cfg, stim_tl_jk_squid);

dippos_squid = stim_dpN20_squid_jk.dip.pos;
           

% dmu3
cfg                 = [];
cfg.sourcemodel.pos = pos275;
cfg.grad            = ctf275;
cfg.headmodel       = headmodel_sphere;
leadfield_squid     = ft_prepare_leadfield(cfg);

sel275 = startsWith(ctf275.label, 'M');
leadfield_squid.leadfield{1,1}=leadfield_squid.leadfield{1,1}(sel275,:);

dipmom_squid = pinv(leadfield_squid.leadfield{1,1}*orient)*stim_dpN20_squid_jk.Vdata ;
   
dmu_squid = sqrt(sum((dipmom_squid-norm_mom_squid).^2)/19);


%% Plot dmu - number of sensors
%% log-log
figure;
plot(log10(nr_sens), log10(dmu/dmu_squid), 'o','MarkerFaceColor', 'RED');
hold on
yline(0, '--')
xlabel('$\log_{10}(N)$','Interpreter', 'latex', 'FontSize', 14);
ylabel('$log_{10}(\sigma_{\hat{q}, OPM}$ / $\sigma_{\hat{q}, SQUID}$)', 'Interpreter', 'latex', 'FontSize', 14);
title('Experiment')

%%% fit a line to the log-log plot.
coefficients = polyfit(log10(nr_sens), log10(dmu/dmu_squid), 1); % Linear fit (degree 1)
a = coefficients(1);
b = coefficients(2);

tmp = linspace(30,510,10);
fitted_line=a*log10(tmp)+b;
plot(log10(tmp), fitted_line, '--');

fitted_line_eq1 = sprintf('Fitted line: y = %.2fx + %.2f', a, b);
legend(fitted_line_eq1, 'Location', 'northwest');

