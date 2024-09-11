addpath('H:/common/matlab/fieldtrip')

ft_defaults

%% ft_databrowser

cfg                = [];
cfg.dataset        = '20240305_143835_sub-001_file-MedianNerve06_raw.fif';
cfg.ylim           = [-1 1] *1e-11; % For EOG/EEG channels cfg.ylim=[-50e-6 50e-6] (50 micro Volt) is a useful scale and for the MEG channels cfg.ylim=[-1e-12 1e-12] (1 pT).
cfg.continuous     = 'yes';
cfg.preproc.demean = 'yes'; % if this is 'no' then the scaling at the buttefly plot is not good. Note that the baselinewindow is automatically selected as the complete 1 sec trial (not the whole 250 sec recording)
cfg.viewmode       = 'butterfly';
cfg.layout         = 'fieldlinebeta2bz_helmet'; % used for the topoplotER
cfg.channel        = {'all', '-ai61'};
ft_databrowser(cfg)

%% Preprocessing for 1 run

event = ft_read_event('20240305_142537_sub-001_file-MedianNerve05_raw.fif');

cfg                    = [];
cfg.dataset            = '20240305_142537_sub-001_file-MedianNerve05_raw.fif';
cfg.trialdef.eventtype = 'ai61';
cfg.trialdef.prestim   = 0.2; % ISI varied from 0.8 to 1.2 sec
cfg.trialdef.poststim  = 0.4;
cfg                    = ft_definetrial(cfg);

cfg.demean         = 'yes';  
cfg.baselinewindow = [-Inf 0]; % see https://www.fieldtriptoolbox.org/workshop/oslo2019/introduction/. cfg.baselinewindow = [-Inf 0] is the same as: [start of the trial, 0] = [-prestim, 0]
cfg.channel        = {'all', '-ai61'};
cfg.continuous     = 'yes'; % without that I get the error: 'requested data segment extends over a discontinuous trial boundary'. I guess I need cfg.continuous since 'subj001ses002_3031000.01_20231208_01.ds' is initially a continuous recording and I do not want to lose the ITI during preprocessing.
data_meg           = ft_preprocessing(cfg);

data_meg.sampleinfo = full(data_meg.sampleinfo);
data_meg.trialinfo  = full(data_meg.trialinfo);

%% Preprocessing for all 6 runs

f = dir('/home/megmethods/kontsi/Documents/recordings/opm-MN/sub-001/*.fif');
f = struct2cell(f);

 for i=1:6
    
    cfg                    = [];
    cfg.dataset            = f{1,i};
    cfg.trialdef.eventtype = 'ai61';
    cfg.trialdef.prestim   = 0.2; % ISI varied from 0.8 to 1.2 sec
    cfg.trialdef.poststim  = 0.4;
    cfg                    = ft_definetrial(cfg);
    
    cfg.demean         = 'yes';  
    cfg.baselinewindow = [-Inf 0]; % see https://www.fieldtriptoolbox.org/workshop/oslo2019/introduction/. cfg.baselinewindow = [-Inf 0] is the same as: [start of the trial, 0] = [-prestim, 0]
    cfg.channel        = {'all', '-ai61'};
    cfg.continuous     = 'yes'; % without that I get the error: 'requested data segment extends over a discontinuous trial boundary'. I guess I need cfg.continuous since 'subj001ses002_3031000.01_20231208_01.ds' is initially a continuous recording and I do not want to lose the ITI during preprocessing.
    data_meg(i)        = ft_preprocessing(cfg);

    data_meg(i).sampleinfo = full(data_meg(i).sampleinfo);
    data_meg(i).trialinfo  = full(data_meg(i).trialinfo);
    
 end

%% HFC

% clear the duplicate channels
data_meg(2).label(find(strcmp(data_meg(2).label, {'R105_bz-1'}))) = {'temp'};
cfg = [];
cfg.channel = {'all' '-R105_bz','-temp'};
data_meg(2) = ft_selectdata(cfg, data_meg(2));

data_meg(4).label(find(strcmp(data_meg(4).label, {'R603_bz-1'}))) = {'temp'};
cfg = [];
cfg.channel = {'all' '-R603_bz','-temp'};
data_meg(4) = ft_selectdata(cfg, data_meg(4));


% rename R101 in run 6
% WAY 1: not working
montage_pos1 = [];
montage_pos1.labelold = data_meg(6).label;

data_meg(6).label(1) = {'R101_bz'};
data_meg(6).label(6) = {'L101_bz'};

montage_pos1.labelnew = data_meg(6).label;
montage_pos1.tra = eye(31);

cfg = [];
cfg.montage = montage_pos1;
data_meg(6) = ft_preprocessing(cfg, data_meg(6));

% WAY 2
data_meg(6).label(find(strcmp(data_meg(6).label, {'R101_bz-1'}))) = {'L101_bz'};

data_meg(6).grad.label{6, 1} = 'L101_bz';


% do HFC
for i=1:6
    data_meg_nohdr(i) = rmfield(data_meg(i), 'hdr'); % Only that way cfg.residualcheck ='yes' works. Because ft_denoise_hfc() -> ft_chanunit() does not include 'ai61' and finds 'Tesla' as a unit. 

    cfg         = [];
    cfg.residualcheck = 'yes';
    data_meg_hfc(i) = ft_denoise_hfc(cfg, data_meg_nohdr(i));
end

%% rejectvisual and timelockanalysis

for i=1:6
    cfg        = [];
    cfg.method = 'summary';
    data_meg_hfc_clean(i) = ft_rejectvisual(cfg, data_meg_hfc(i)); % ft_rejectvisual: reject all trials that have more than 5*10^-23 T^2 variance VS ft_rejectvisual + HFC: reject all trials that "stick out" (~ 5*10^-23 T^2 variance) 
    
    cfg     = [];
    stim_tl(i) = ft_timelockanalysis(cfg, data_meg_hfc_clean(i));
end

% cfg = [];
% cfg.layout = 'fieldlinebeta2bz_helmet';
% cfg.xlim = [0 0.2];
% cfg.zlim = [-1 1] *1e-12;
% ft_movieplotER(cfg, stim_tl(1)); 

%% ft_topoplotER for all sensors in each run
for i=1:6
    cfg = [];
    cfg.layout = 'fieldlinebeta2bz_helmet.mat';
    cfg.skipcomnt = 'no';
    cfg.skipscale = 'no';
    layout = ft_prepare_layout(cfg);

    cfg = [];
    cfg.layout = rmfield(layout, 'mask');
    cfg.channel = stim_tl(i).label;
    cfg.mask = 'convex';
    cfg.outline = 'convex'; % I need to specify that. Otherwise the outline is too big in comparison to the sensors locations
    layout_trimmed = ft_prepare_layout(cfg);

    cfg = [];
    cfg.xlim = [0.035 0.050];
    cfg.layout = layout_trimmed;
    cfg.zlim = 'maxabs';
    ft_topoplotER(cfg, stim_tl(i)); 

    filename = sprintf('run%d.jpeg', i);
    saveas(gcf, filename);
end 

close all

%% Power spectral density (PSD)
% PSD: demean & ft_rejectvisual & HFC
cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [0 600];
cfg.keeptrials = 'no'; % I plot the average PSD for every sensor averaged across all trials
cfg.tapsmofrq  = 2;
freq       = ft_freqanalysis(cfg, data_meg_hfc_clean(1));

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
freq       = ft_freqanalysis(cfg, data_meg_clean(1));

figure; plot(freq.freq,  (log10(freq.powspctrm))); hold on; plot(freq.freq, mean(log10(freq.powspctrm)), 'k', 'LineWidth', 2)
legend(freq.grad.label)

%% Check how much the subject moved its head by comparing the 9 stable sensors: not recommended

%% Is the topo the same for each sensor in each run (if not the stimulator's timing could have changed between runs)

% use the 9 stable sensors


% A. choose a time point of interest (35 -50 ms) and compare the
% topographies for these 9 sensors only for every experimental run (so 6
% topoplotER in total). Note to avoid interpolation use a trimmed mask as
% in: https://www.fieldtriptoolbox.org/tutorial/preprocessing_opm/
for i = 1:6
    cfg = [];
    cfg.layout = 'fieldlinebeta2bz_helmet.mat';
    cfg.skipcomnt = 'no';
    cfg.skipscale = 'no';
    layout = ft_prepare_layout(cfg);

    cfg = [];
    cfg.layout = rmfield(layout, 'mask');
    cfg.channel = {'R101_bz', 'L108_bz', 'L503_bz', 'R503_bz', 'L507_bz', 'R507_bz', 'L212_bz', 'R212_bz'};
    cfg.mask = 'convex';
    cfg.outline = 'convex'; % I need to specify that. Otherwise the outline is too big in comparison to the sensors locations
    layout_trimmed = ft_prepare_layout(cfg);

    cfg = [];
    cfg.xlim = [0.035 0.050];
    cfg.layout = layout_trimmed;
    cfg.marker = 'labels';
    cfg.zlim = 'maxabs';
    ft_topoplotER(cfg, stim_tl(i)); 
end

% B. plot the singleplotER for 1 sensor * 6 runs (so 1*6 singleplotER in
% total)
for i = 1:6
    cfg = [];
    cfg.channel = {'R507_bz'};
    cfg.xlim = [-0.1 0.2];
    cfg.ylim = [-5e-13, 5e-13];
    ft_singleplotER(cfg, stim_tl(i));

%     % Save the figure in the current folder
%     filename = sprintf('run%d.jpeg', i);
%     saveas(gcf, filename);
%     close(figure); % Close the figure window to prevent accumulation
end




%% Are the sensors in the right slots?
%%

% remove -
for k=1:6
    for i = 1:numel(stim_tl(k).label)
      tok = tokenize(stim_tl(k).label{i}, '-');
      stim_tl(k).label{i} = tok{1};
    end
end

%%

% find channels that might appear twice due to wrong labelling - this
% problem occured during opm installation. Then delete these these
% duplicate channels. 
% 
% NOTE THIS STEP IS NEEDED FOR "only demeaning" AND
% "ft_rejectvisual" BUT NOT FOR "hfc". For "hfc" I HAVE AREADY DONE IT IN
% data_meg_hfc

% Do that only for the first 5 runs
for k=1:6
    badchannel = [];
    for i=1:length(stim_tl(k).label)
      for j=(i+1):length(stim_tl(k).label)
        if strcmp(stim_tl(k).label{i}, stim_tl(k).label{j})
          badchannel(end+1) = i;
          badchannel(end+1) = j;
        end
      end
    end

    % Print channels with the same name
    fprintf('Channels with the same name in run %d:\n', k);
    for b = 1:length(badchannel)
        fprintf('%s\n', stim_tl(k).label{badchannel(b)});
    end
    fprintf('\n');
end

% Channels with the same name in run 1:
% 
% Channels with the same name in run 2:
% R105_bz
% R105_bz
% 
% Channels with the same name in run 3:
% 
% Channels with the same name in run 4:
% R603_bz
% R603_bz
% 
% Channels with the same name in run 5:
% 
% Channels with the same name in run 6:
% R101_bz
% R101_bz

% clear the duplicate channels in run 2 and 4
stim_tl(2).label(19) = {'temp'};
cfg = [];
cfg.channel = {'all' '-R105_bz','-temp'};
stim_tl(2) = ft_timelockanalysis(cfg, stim_tl(2));

stim_tl(4).label(15) = {'temp'};
cfg = [];
cfg.channel = {'all' '-R603_bz','-temp'};
stim_tl(4) = ft_timelockanalysis(cfg, stim_tl(4));

% In run 6, R101_bz appears also twice. However, since this is one of
% the stable sensors I can check where it was placed in the 5th run. See
% header5.orig.ch_names:

% header5 = ft_read_header('20240305_142537_sub-001_file-MedianNerve05_raw.fif');
% header6 = ft_read_header('20240305_143835_sub-001_file-MedianNerve06_raw.fif');

% So, in the R101 slot is the s1 sensor. Now I look at my notes and I see
% that the other R101_bz-s6 sensor should go to L101, since I expect this
% slot in run6 but it is not there
stim_tl(6).label(1) = {'R101_bz'};
stim_tl(6).label(6) = {'L101_bz'};

stim_tl(6).grad.label{6, 1} = 'L101_bz';

%%

% I had 8 stable sensors between runs. I will rename them
for k = 2:6
    index = find(ismember(stim_tl(k).label, {'R101_bz', 'L108_bz','R503_bz','L503_bz','L507_bz','R507_bz','R212_bz','L212_bz'}));
    for i = 1:length(index)
        stim_tl(k).label{index(i)} = strcat(stim_tl(k).label{index(i)}, '_', num2str(k));
    end
end

%%

% the 6th run has some channels that appeared in the 5th run too. I am going to
% rename it. 
index = find(ismember(stim_tl(6).label, stim_tl(5).label));

for i = 1:length(index)   
        stim_tl(6).label{index(i)} = strcat(stim_tl(6).label{index(i)}, '-', num2str(6));
end

% Also R114 appears in the 6th and 1st run so I rename it.
index_R114 = find(strcmp(stim_tl(6).label, 'R114_bz'));

if ~isempty(index_R114)
    stim_tl(6).label{index_R114} = 'R114_bz-6';
end

%%

% After I checked all the channels that I know that they are duplicate, now I will see which channels still remain duplicate
duplicate_channels = {};

for i = 1:6
    % Loop through each channel in the current run
    for j = 1:length(stim_tl(i).label)
        current_channel = stim_tl(i).label{j};
        % Compare the current channel with channels in other runs
        for k = (i + 1):6
            % Check if the channel name exists in the other runs
            if any(strcmp(current_channel, stim_tl(k).label))
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

% Duplicate Channels:
% Channel Name       :   Runs
% L214_bz   :   Run 2 and Run 3
% R504_bz   :   Run 3 and Run 5
% L504_bz   :   Run 3 and Run 5

% Now I rename these duplicate channels with the logic: L214_bz is
% duplicate at Run 2 and Run 3. So, I rename the run 3 as L214_bz__3:
stim_tl(3).label{27}=strcat(stim_tl(3).label{27}, '__3');
stim_tl(5).label{16}=strcat(stim_tl(5).label{16}, '__5');
stim_tl(5).label{31}=strcat(stim_tl(5).label{31}, '__5');


% CONCLUSION

% there is no L105, two R105 appear in the 2nd run (problem with the
% HEDscan labeling)

% two R603 in the 4th run

% L312 missing (probably human mistake during the experiment. Mixed it with L214)

% R/L604 probably missing too (probably human mistake during the
% experiment. Mixed it with R/L504)

% So in the end I expect 7 out of 144 channels missing.


%% Append

cfg       = [];
cfg.appenddim = 'chan';
append_tl = ft_appendtimelock(cfg, stim_tl(1), stim_tl(2), stim_tl(3), stim_tl(4), stim_tl(5), stim_tl(6));
% Warning: renaming the appended averages to "trial" 

cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.showlabels = 'yes';
cfg.ylim = 'maxabs';
ft_multiplotER(cfg, append_tl);

cfg = [];
cfg.xlim = [0.035 0.050];
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.zlim = 'maxabs';
ft_topoplotER(cfg, append_tl); 

% R101_bz, R101_bz_2, R101_bz_3 ... are not averaged by ft_multiplotER.
% ft_multiplotER keeps only the ERF of R101_bz. I know that since if I
% delete R101_bz from append_tl then no channel appears at the location of
% R101_bz when I use ft_multiplotER.

% from ft_multiplotER I see that I miss 7 channels in the end:
% R/L105, L312, R/L604, R603, R605

 