
fl=ft_read_sens('fieldlinebeta2.mat')
fl=ft_convert_units(fl, 'mm')


% headmodel = ft_read_vol('standard_singleshell.mat');  % Load a standard BEM model
% headmodel=ft_convert_units(headmodel, 'mm')

singlesphere = [];
singlesphere.r = 100; % singlesphere.r = 0.1*10^2  
singlesphere.o = [0 0 0]; % singlesphere.o = [0.01*10^2 0 0] 
singlesphere.unit = 'mm';

cfg.dip.pos = [0 0 60]; % left motor cortex
cfg.dip.mom = [0 -1 0]; 
figure
ft_plot_headmodel(singlesphere, 'facealpha', 0.1, 'facecolor', 'brain')
ft_plot_sens(fl)
ft_plot_dipole(cfg.dip.pos,cfg.dip.mom, 'unit', 'mm')
view([90 0]); % Lateral view
camlight
lighting gouraud
print('1', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution

cfg.dip.pos = [0 0 60]; % left motor cortex
cfg.dip.mom = [0 -1 0]; 
figure
ft_plot_headmodel(singlesphere, 'facealpha', 0.1, 'facecolor', 'brain')
ft_plot_sens(fl)
ft_plot_dipole(cfg.dip.pos,cfg.dip.mom, 'unit', 'mm')
view([90 90]); % Lateral view
camlight
lighting gouraud
print('1', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution

cfg = [];
cfg.ntrials = 10;
cfg.triallength = 10;
cfg.fsample = 1000;
cfg.dip.unit = 'mm';
cfg.dip.pos = [0 0 60]; 
cfg.dip.mom = [0 -1 0]; 
cfg.frequency = 0;
cfg.grad = fl;
cfg.headmodel = singlesphere;
data = ft_dipolesimulation(cfg);

tl=ft_timelockanalysis([],data);

% Invert the data values by multiplying by -1
tl.avg = tl.avg * -1;  

cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.channel         = {'all'}; 
cfg.markersize = 3;
cfg.colormap = flipud(brewermap([], 'RdBu'));  % Invert the 'RdBu' colormap
% cfg.colormap= 'RdBu';
ft_topoplotER(cfg, tl)
colorbar

print('2', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution
