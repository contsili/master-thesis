 cfg                    = [];
    cfg.dataset            = 'M:\Documents\recordings\squid-MN\empty-room\emptyroom_Noise_20231207_01.ds';
%     cfg.bsfilter = 'yes';
%     cfg.bsfreq=[49 51; 99 101; 149 151; 199 201; 249 251; 299 301];
    cfg.demean             = 'yes';  
    cfg.baselinewindow     = 'all'; % see https://www.fieldtriptoolbox.org/workshop/oslo2019/introduction/. cfg.baselinewindow = [-Inf 0] is the same as: [start of the trial, 0] = [-prestim, 0]
    cfg.channel            = {'all'};
    cfg.continuous         = 'yes'; % without that I get the error: 'requested data segment extends over a discontinuous trial boundary'. I guess I need cfg.continuous since 'subj001ses002_3031000.01_20231208_01.ds' is initially a continuous recording and I do not want to lose the ITI during preprocessing.
    data_meg            = ft_preprocessing(cfg);

cfg = [];
    cfg.length  = 1;
    cfg.overlap = 0;
    data_meg = ft_redefinetrial(cfg, data_meg);

    sel275 = startsWith(data_meg.label, 'M');

    cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [10 300];
cfg.keeptrials = 'no'; % average across trials
cfg.tapsmofrq  = 2;
% cfg.channel = data_meg.label(sel275,:);
cfg.channel =  {'MLC21', 'MLC32', 'MLC41', 'MLC51', 'MLC52', 'MLC53', 'MLC54', 'MLC55', 'MLC61', 'MLC62', 'MLC63', ...
                       'MLO11', 'MLO12', 'MLO13', 'MLO21', 'MLO22', 'MLP11', 'MLP12', 'MLP21', 'MLP22', 'MLP23', 'MLP31', ...
                       'MLP32', 'MLP34', 'MLP35', 'MLP41', 'MLP42', 'MLP43', 'MLP51', 'MLP52', 'MLP54', 'MRC13', 'MRC14', ...
                       'MRC15', 'MRC16', 'MRC17', 'MRC21', 'MRC22', 'MRC23', 'MRC24', 'MRC25', 'MRC31', 'MRC32', 'MRC41', ...
                       'MRC42', 'MRC51', 'MRC52', 'MRC53', 'MRC54', 'MRC55', 'MRC61', 'MRC63', 'MRF46', 'MRF55', 'MRF56', ...
                       'MRF64', 'MRF65', 'MRF66', 'MRF67', 'MRO11', 'MRO12', 'MRO13', 'MRO14', 'MRO21', 'MRO22', 'MRO23', ...
                       'MRO24', 'MRO34', 'MRP11', 'MRP12', 'MRP21', 'MRP22', 'MRP23', 'MRP31', 'MRP32', 'MRP33', 'MRP34', ...
                       'MRP35', 'MRP41', 'MRP42', 'MRP43', 'MRP44', 'MRP45', 'MRP51', 'MRP52', 'MRP53', 'MRP54', 'MRP55', ...
                       'MRP56', 'MRP57', 'MRT13', 'MRT14', 'MRT15', 'MRT16', 'MRT25', 'MRT26', 'MRT27', 'MRT37', 'MZC02', ...
                       'MZC03', 'MZC04', 'MZO01', 'MZP01'};
freq       = ft_freqanalysis(cfg, data_meg);

freq.powspctrm_ft_sqrtHz = sqrt(freq.powspctrm) * 1e15;

hold on;
loglog(freq.freq, freq.powspctrm_ft_sqrtHz, 'LineWidth', 0.001, 'Color', [.5 0 .5]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (fT/\surdHz)');
grid on;

set(gca, 'FontSize', 12);

ax = gca;
ax.YTick = [5,10,50, 100]; % Define y-axis tick positions

ylim([2,100])
print('e3', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution
