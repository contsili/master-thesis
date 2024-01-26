%% Magnetic dipole phantom



cfg = [];
cfg.dataset = 'C:\Users\user\Documents\Courses\Internship\master-thesis\phantom\raw\01_3031000.01_20231010_01.ds';
cfg.preproc.demean = 'yes';
cfg.ylim = [-1 1]*1e-11;
ft_databrowser(cfg)


%% 

% hdr = ft_read_header('E:\3031000.01\raw\01_3031000.01_20231010_01.ds');
% 
% event = ft_read_event('E:\3031000.01\raw\01_3031000.01_20231010_01.ds');


%%

cfg = [];
cfg.dataset = 'C:\Users\user\Documents\Courses\Internship\master-thesis\phantom\raw\01_3031000.01_20231010_01.ds';
cfg.detrend = 'yes';

cfg.lpfilter  = 'yes';
cfg.lpfreq    = 40;
cfg.lpfiltord = 2;
data = ft_preprocessing(cfg);


%%
cfg = [];
timelock = ft_timelockanalysis(cfg, data); % sampling rate is 600 and not 625 as it says in the CTF guide

%%
% timelock.avg(1,:) = []; % ignore SCLK01 channel
% figure;
% plot(timelock.time, mean(timelock.avg));
% title('Average across 368 CTF channels')

%%
hold off

cfg = [];
cfg.layout = 'C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\layout\CTF275_helmet';
cfg.colormap = '-RdBu';
ft_multiplotER(cfg, timelock);

%%
figure;
cfg = [];
cfg.layout = 'C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\layout\CTF275_helmet';
cfg.colormap = '-RdBu';
ft_movieplotER(cfg, timelock)

%%
% cfg = [];
% cfg.method = 'infinite';
% headmodel = ft_prepare_headmodel(cfg); 
% 
% ft_plot_headmodel(headmodel) % Warning: there is nothing to plot for an infinite volume conductor 


%%
 load 'C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\ctf275.mat';
% % neuromag306 = ft_convert_units(neuromag306, 'cm');
% figure;
% ft_plot_sens(ctf275, 'label', 'no', 'axes', 1, 'orientation', 1, 'marker','o');

%%

cfg=[];
cfg.vol      = [];
cfg.vol.type = 'infinite';
cfg.numdipoles=1;
cfg.model       = 'moving';
cfg.gridsearch = 'no';
cfg.nonlinear = 'yes';
cfg.dip.pos = [0,0,0.05];
cfg.grad             = ctf275;
% cfg.senstype        = 'meg';            % sensor type
source = ft_dipolefitting(cfg, timelock);

%%

figure;

for i=1:length(source.dip)

    disp(i);
    cla
    
    ft_plot_sens(ctf275, 'label', 'no', 'axes', 0, 'orientation', 0, 'marker','o', 'unit', 'cm'); % I have to specify ('unit', 'cm') in order for the dipole and the sensors to have the same size 
    hold on
    ft_plot_dipole(source.dip(i).pos,source.dip(i).mom);

%     view(120,0); 

    drawnow 
end

% It is beautiful. It plots the magnetic dipole (not the electric one!)


%% TFR: See the 7 Hz component. Also what else can I do? 1. I can play with my 2nd recording where I have 100 trials (not 1)
