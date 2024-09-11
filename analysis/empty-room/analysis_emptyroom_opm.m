addpath('H:/common/matlab/fieldtrip')

ft_defaults

%% Preprocessing for all 6 runs

f = dir('M:\Documents\recordings\opm-MN\empty-room\*.fif');
f = struct2cell(f);

 for i=1:6
    
    cfg                    = [];
    cfg.dataset            = f{1,i};
    cfg.demean             = 'yes';  
    cfg.baselinewindow     = 'all'; % see https://www.fieldtriptoolbox.org/workshop/oslo2019/introduction/. cfg.baselinewindow = [-Inf 0] is the same as: [start of the trial, 0] = [-prestim, 0]
    cfg.channel            = {'all', '-ai61'};
    cfg.continuous         = 'yes'; % without that I get the error: 'requested data segment extends over a discontinuous trial boundary'. I guess I need cfg.continuous since 'subj001ses002_3031000.01_20231208_01.ds' is initially a continuous recording and I do not want to lose the ITI during preprocessing.
    data_meg(i)            = ft_preprocessing(cfg);

 end


 % In the emptyroom recording all the sensors were working. The additional sensor
 % was put in L201 for run 6 and L101 for run 1-5


% create 1 sec "fake" trials. Otherwise PSD calculation takes too long 
for i=1:6
    cfg = [];
    cfg.length  = 1;
    cfg.overlap = 0;
    data_meg(i) = ft_redefinetrial(cfg, data_meg(i));
end


%% HFC

% No duplicate channels
% do HFC
for i=1:6
    data_meg_nohdr(i) = rmfield(data_meg(i), 'hdr'); % Only that way cfg.residualcheck ='yes' works. Because ft_denoise_hfc() -> ft_chanunit() does not include 'ai61' and finds 'Tesla' as a unit. 

    cfg         = [];
    cfg.residualcheck = 'yes';
    data_meg_hfc(i) = ft_denoise_hfc(cfg, data_meg_nohdr(i));
end

%% I do not do ft_rejectvisual since I have not split my data in trials

%% Power spectral density (PSD)

% NOTE data_meg_hfc(1) corresponds to run6!!

% PSD: demean & ft_rejectvisual & HFC
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [0 100];
cfg.keeptrials = 'no'; % I plot the average PSD for every sensor averaged across all trials
cfg.tapsmofrq  = 2;
freq       = ft_freqanalysis(cfg, data_meg_hfc(1));

figure; plot(freq.freq,  (log10(freq.powspctrm))); hold on; plot(freq.freq, mean(log10(freq.powspctrm)), 'k', 'LineWidth', 2)
legend(freq.grad.label)

% A channel remeains noisy for high freq. (>300 Hz). Find which channel it
% is:
g = log10(freq.powspctrm(:,freq.freq==400)); 
index = find(g==max(g)); 
freq.grad.label(index); % it is R503



% PSD: demean & ft_rejectvisual
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [0 100];
cfg.keeptrials = 'no'; % I plot the average PSD for every sensor averaged across all trials
cfg.tapsmofrq  = 2;
freq       = ft_freqanalysis(cfg, data_meg(1));

figure; plot(freq.freq,  (log10(freq.powspctrm))); hold on; plot(freq.freq, mean(log10(freq.powspctrm)), 'k', 'LineWidth', 2)
legend(freq.grad.label)





% PSD of append_data_resample
load M:\Documents\recordings\opm-MN\empty-room\preprocessing\append_data_resample.mat
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [10 300];
cfg.keeptrials = 'no'; % average across trials
cfg.tapsmofrq  = 2;
cfg.channel = {'R304_bz','R403_bz','R210_bz', 'R306_bz', 'R503_bz', 'R207_bz', 'R408_bz', 'R305_bz', 'R209_bz', 'R205_bz', 'R407_bz', 'R505_bz', 'R208_bz', 'R309_bz', 'R406_bz', 'R206_bz', 'R308_bz', 'R404_bz', 'R307_bz'};

freq       = ft_freqanalysis(cfg, append_data_resample);

freq.powspctrm_ft_sqrtHz = sqrt(freq.powspctrm) * 1e15;

hold on;
loglog(freq.freq, freq.powspctrm_ft_sqrtHz, 'LineWidth', 0.5, 'Color', [.5 0 .5]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (fT/\surdHz)');
grid on;

ax = gca;
ax.YTick = [5,10,50, 100]; % Define y-axis tick positions

set(gca, 'FontSize', 12);
ylim([2,100])

print('e1', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution


%% Noise covariance

% no hfc
cfg=[];
cfg.resamplefs = 600;
data_meg_resample(1) = ft_resampledata(cfg, data_meg(1));

cfg                  = [];
cfg.covariance       = 'yes';
timelock_nohfc       = ft_timelockanalysis(cfg, data_meg_resample(1));


% hfc
cfg=[];
cfg.resamplefs = 600;
data_meg_hfc_resample(1) = ft_resampledata(cfg,  data_meg_hfc(1));

cfg                  = [];
cfg.covariance       = 'yes';
timelock_hfc         = ft_timelockanalysis(cfg, data_meg_hfc_resample(1));


% compare off-diagonal elements
off_diagonal = timelock_nohfc.cov(~eye(size(timelock_nohfc.cov))); % Extract all non-diagonal elements
mean_off_diagonal_nohfc = mean(off_diagonal) % Compute the mean

off_diagonal = timelock_hfc.cov(~eye(size(timelock_hfc.cov))); % Extract all non-diagonal elements
mean_off_diagonal_hfc = mean(off_diagonal) % Compute the mean

% mean_off_diagonal_nohfc =
% 
%    1.1815e-24
% 
% 
% mean_off_diagonal_hfc =
% 
%    2.0309e-25





%% Noise covariance - diagonal elements

%% Rename sensor that re-appear between runs -
% 
% !!!!! Remember i=1 for the empty room corresponds to i=6 for the MN experiment !!!!!

% I had 8 stable sensors between runs. I will rename them
for i = 1:5     
    index = find(ismember(data_meg_hfc(i).label, {'R101_bz', 'L108_bz','R503_bz','L503_bz','L507_bz','R507_bz','R212_bz','L212_bz'}));
    for k = 1:length(index)
        data_meg_hfc(i).label{index(k)} = strcat(data_meg_hfc(i).label{index(k)}, '_', num2str(i));
    end    
end


% rename L101_bz
for i=1:5
index = find(ismember(data_meg_hfc(i).label, {'L101_bz'}));

    for k = 1:length(index)
        data_meg_hfc(i).label{index(k)} = strcat(data_meg_hfc(i).label{index(k)}, '_', num2str(i));
    end 
end


% rename 'R114_bz', 'L201_bz'
for i=1:5
index = find(ismember(data_meg_hfc(i).label, {'R114_bz', 'L201_bz'}));

    for k = 1:length(index)
        data_meg_hfc(i).label{index(k)} = strcat(data_meg_hfc(i).label{index(k)}, '_', num2str(i));
    end 
end

% rename 'L504_bz'
for i=1:2
index = find(ismember(data_meg_hfc(i).label, {'L504_bz'}));

    for k = 1:length(index)
        data_meg_hfc(i).label{index(k)} = strcat(data_meg_hfc(i).label{index(k)}, '_', num2str(i));
    end 
end

% rename R504 from run 4 (which is in the position i=2)
for i=2
index = find(ismember(data_meg_hfc(i).label, {'R504_bz'}));

    for k = 1:length(index)
        data_meg_hfc(i).label{index(k)} = strcat(data_meg_hfc(i).label{index(k)}, '_', num2str(i));
    end 
end


%% After I checked all the channels that I know that they are duplicate, now I will see which channels still remain duplicate
duplicate_channels = {};

for i = 1:6
    % Loop through each channel in the current run
    for j = 1:length(data_meg_hfc(i).label)
        current_channel = data_meg_hfc(i).label{j};
        % Compare the current channel with channels in other runs
        for k = (i + 1):6
            % Check if the channel name exists in the other runs
            if any(strcmp(current_channel, data_meg_hfc(k).label))
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

%%
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\squid_noise_floor.mat


cfg         = [];
cfg.keepsampleinfo='no';
append_data = ft_appenddata(cfg, data_meg_hfc(1), data_meg_hfc(2), data_meg_hfc(3), data_meg_hfc(4), data_meg_hfc(5), data_meg_hfc(6)); 

% resample cause the file is too big and ft_timelockanalysis does not work
cfg=[];
cfg.resamplefs = 1200;
append_data_resample = ft_resampledata(cfg,  append_data);

%%
cfg             = [];
cfg.covariance       = 'yes';
append_tl    = ft_timelockanalysis(cfg, append_data_resample);

diagonal = diag(append_tl.cov);
diagonal_std = sqrt(diagonal);

append_tl.ratio = diagonal_std/squid_noise_floor;

figure
cfg = [];
cfg.xlim = [0,0];
cfg.parameter = 'ratio';
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.marker = 'labels';
cfg.channel={'all','-L101_bz','-L105_bz', '-R105_bz', '-L214_bz','-L312_bz', '-R603_bz', '-R604_bz', '-R605_bz', '-L604_bz'}
ft_topoplotER(cfg, append_tl); 
colorbar
% hcb = colorbar;
% hcb.Ruler.Exponent = -13;

mean(append_tl.ratio)

% mean(append_tl.ratio)
% 
% ans =
% 
%     3.6734

% Note sensor R505 has 2.4 more noise than SQUID.


%% mean noise ratio on the dipolar pattern

% yellow: R210_bz, R209_bz, R208_bz, R309_bz, R308_bz, R307_bz
% blue: R306_bz, R405_bz, R504_bz, R305_bz, R505_bz, R406_bz, R404_bz

% try2: R210_bz, R306_bz, R503_bz, R405_bz, R207_bz, R504_bz, R408_bz, R305_bz, R209_bz, R205_bz, R407_bz, R505_bz, R403_bz, R208_bz, R309_bz, R406_bz, R304_bz, R206_bz, R308_bz, R404_bz, R307_bz

% Assuming append_tl_jk(1).label is a cell array of labels
labels = append_tl.label;

% Names to search for (ignore R504 and R405)
names_to_find = {'R304_bz','R403_bz','R210_bz', 'R306_bz', 'R503_bz', 'R207_bz', 'R408_bz', 'R305_bz', 'R209_bz', 'R205_bz', 'R407_bz', 'R505_bz', 'R208_bz', 'R309_bz', 'R406_bz', 'R206_bz', 'R308_bz', 'R404_bz', 'R307_bz'};

% Initialize an empty cell array to store indices of matching names
matching_indices = {};

% Iterate through labels and find indices of matching names
for idx = 1:numel(labels)
    if any(strcmp(labels{idx}, names_to_find))
        matching_indices{end+1} = idx;
    end
end

% Convert the cell array to a numeric array if needed
matching_indices = cell2mat(matching_indices);

opm_noise_floor_opm_dippat = mean(diagonal_std(matching_indices,:))

% OPM vs the whole SQUID array
mean(append_tl.ratio(matching_indices,:))

% ans =
% 
%     3.2134

% OPM vs the SQUID array on top of the dipolar pattern
names_to_find_squid = {'MLC21', 'MLC32', 'MLC41', 'MLC51', 'MLC52', 'MLC53', 'MLC54', 'MLC55', 'MLC61', 'MLC62', 'MLC63', ...
                       'MLO11', 'MLO12', 'MLO13', 'MLO21', 'MLO22', 'MLP11', 'MLP12', 'MLP21', 'MLP22', 'MLP23', 'MLP31', ...
                       'MLP32', 'MLP34', 'MLP35', 'MLP41', 'MLP42', 'MLP43', 'MLP51', 'MLP52', 'MLP54', 'MRC13', 'MRC14', ...
                       'MRC15', 'MRC16', 'MRC17', 'MRC21', 'MRC22', 'MRC23', 'MRC24', 'MRC25', 'MRC31', 'MRC32', 'MRC41', ...
                       'MRC42', 'MRC51', 'MRC52', 'MRC53', 'MRC54', 'MRC55', 'MRC61', 'MRC63', 'MRF46', 'MRF55', 'MRF56', ...
                       'MRF64', 'MRF65', 'MRF66', 'MRF67', 'MRO11', 'MRO12', 'MRO13', 'MRO14', 'MRO21', 'MRO22', 'MRO23', ...
                       'MRO24', 'MRO34', 'MRP11', 'MRP12', 'MRP21', 'MRP22', 'MRP23', 'MRP31', 'MRP32', 'MRP33', 'MRP34', ...
                       'MRP35', 'MRP41', 'MRP42', 'MRP43', 'MRP44', 'MRP45', 'MRP51', 'MRP52', 'MRP53', 'MRP54', 'MRP55', ...
                       'MRP56', 'MRP57', 'MRT13', 'MRT14', 'MRT15', 'MRT16', 'MRT25', 'MRT26', 'MRT27', 'MRT37', 'MZC02', ...
                       'MZC03', 'MZC04', 'MZO01', 'MZP01'};

labels_squid = stim_tl.label;

% Initialize an empty cell array to store indices of matching names
matching_indices = {};

% Iterate through labels and find indices of matching names
for idx = 1:numel(labels_squid)
    if any(strcmp(labels_squid{idx}, names_to_find_squid))
        matching_indices{end+1} = idx;
    end
end

% Convert the cell array to a numeric array if needed
matching_indices = cell2mat(matching_indices);

squid_noise_floor_dippat = mean(stim_tl.std(matching_indices, 1:240), 'all');

opm_noise_floor_opm_dippat/squid_noise_floor_dippat

% ans =
% 
%     2.8787

% without R304 and R403:
% ans =
% 
%     2.5423
