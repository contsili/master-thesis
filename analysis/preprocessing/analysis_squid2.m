% Robert suggestion: do not filter. Since the onset of my stim trials is
% 20%-jittered the high-freq. RF artifact and the 50 Hz power line will average out when we calculate the ERF. 
% 
% Also here I ignore the omission trials and I fit only the S1 dipole from
% the stim trials
%
% Even after this "less conservative" approach we see that the S1 dipole
% appears at 35-50 ms and not at 20 ms 

cfg                     = [];
cfg.dataset             = 'subj001ses002_3031000.01_20231208_01.ds';
cfg.trialdef.eventtype  = 'UPPT001';
cfg.trialdef.eventvalue = [1 4];
cfg.trialdef.prestim    = 0.2; % ISI varied from 0.8 to 1.2 sec
cfg.trialdef.poststim   = 0.4;
cfg                     = ft_definetrial(cfg);


cfg.demean         = 'yes';  
cfg.baselinewindow = [-Inf 0]; % see https://www.fieldtriptoolbox.org/workshop/oslo2019/introduction/. cfg.baselinewindow = [-Inf 0] is the same as: [start of the trial, 0] = [-prestim, 0]
cfg.channel        = 'MEG';
cfg.continuous     = 'yes'; % without that I get the error: 'requested data segment extends over a discontinuous trial boundary'. I guess I need cfg.continuous since 'subj001ses002_3031000.01_20231208_01.ds' is initially a continuous recording and I do not want to lose the ITI during preprocessing.
data_meg           = ft_preprocessing(cfg);

cfg        = [];
cfg.trials = (data_meg.trialinfo==1);
stim       = ft_selectdata(cfg, data_meg);

%%

cfg        = [];
cfg.method = 'summary';
stim_clean = ft_rejectvisual(cfg,stim);

%%
cfg     = [];
stim_tl = ft_timelockanalysis(cfg, stim_clean);

cfg = [];
cfg.layout = 'CTF275_helmet';
cfg.xlim = [0 0.2];
cfg.zlim = [-1 1] *1e-13;
ft_movieplotER(cfg, stim_tl); 

cfg = [];
cfg.xlim = [0.035 0.050];
cfg.layout = 'CTF275_helmet';
cfg.zlim = 'maxabs';
cfg.marker ='labels';
ft_topoplotER(cfg, stim_tl); 

% at 75 ms we see the S1 dipole remaining while the S2 dipole appears
% contralaterally and ipsilaterally

%% dipole fit
% load mri_segmented

mri_segmented = ft_convert_units(mri_segmented, 'cm');

cfg = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented);

cfg.funparameter = 'skull';
ft_sourceplot(cfg, mri_segmented);

cfg.funparameter = 'scalp';
ft_sourceplot(cfg, mri_segmented);

%% ft_prepare_mesh

cfg                = [];
cfg.tissue         = {'brain'};
mesh_brain          = ft_prepare_mesh(cfg, mri_segmented);

cfg                = [];
cfg.tissue         = {'skull'};
mesh_skull         = ft_prepare_mesh(cfg, mri_segmented);

cfg                = [];
cfg.tissue         = {'scalp'};
mesh_scalp         = ft_prepare_mesh(cfg, mri_segmented);

figure;
ft_plot_mesh(mesh_scalp)

%% singlesphere headmodel 

cfg        = [];
cfg.method = 'singlesphere';
headmodel_sphere  = ft_prepare_headmodel(cfg, mesh_brain);

headmodel_sphere = ft_convert_units(headmodel_sphere, 'cm');

%% sensor array

ctf275 = ft_read_sens('subj001ses002_3031000.01_20231208_01.ds', 'senstype', 'meg');

%% Dipole fit

cfg            = [];
cfg.latency    = [0.035 0.05];
cfg.unit       = 'cm';
cfg.gridsearch = 'yes';
cfg.grad       = ctf275;
cfg.headmodel  = headmodel_sphere;
stim_dpN20_squid     = ft_dipolefitting(cfg, stim_tl);

figure
hold on

ft_plot_dipole(stim_dpN20_squid.dip.pos, mean(stim_dpN20_squid.dip.mom(1:3,:),2), 'color', 'r')

pos = mean(stim_dpN20_squid.dip.pos,1);
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off

view(12, -10)


%% COMPARISON: Plot opm and squid dipole in the same figure

% squid 
% Express opm and squid in the same coordinate system
% mri_segmented_squid = ft_convert_coordsys(mri_segmented_squid, 'neuromag');
ctf275 = ft_convert_coordsys(ctf275, 'neuromag'); % neuromag coordinate system = 'ras'

% I should use the OPM headmodel_sphere and mri_segmented which are in 'neuromag' coordsys 
load M:\Documents\recordings\opm-MN\results\coregistration\hfc\'mat files'\headmodel_sphere.mat
load M:\Documents\recordings\opm-MN\results\coregistration\hfc\'mat files'\mri_segmented.mat
load M:\Documents\recordings\squid-MN\results\'demeaning & ft_rejectvisual'\stim_tl.mat

cfg            = [];
cfg.latency    = [0.035 0.05];
cfg.unit       = 'cm';
cfg.gridsearch = 'yes';
cfg.grad       = ctf275;
cfg.headmodel  = headmodel_sphere;
stim_dpN20_squid     = ft_dipolefitting(cfg, stim_tl);

% opm
append_tl.avg = append_tl.trial;

cfg            = [];
cfg.latency    = [0.035 0.05];
cfg.unit       = 'cm';
cfg.gridsearch = 'yes';
cfg.grad       = fieldlinebeta2_head;
cfg.channel    = {'all', '-R105_bz', '-R603_bz', '-L101_bz', '-L214_bz'}; 
cfg.headmodel  = headmodel_sphere;
stim_dpN20_opm     = ft_dipolefitting(cfg, append_tl);

figure
hold on

ft_plot_dipole(stim_dpN20_opm.dip.pos, mean(stim_dpN20_opm.dip.mom(1:3,:),2), 'color', 'b')
ft_plot_dipole(stim_dpN20_squid.dip.pos, mean(stim_dpN20_squid.dip.mom(1:3,:),2), 'color', 'r')

legend

% plot the opm crosshair. So load mri_segmented by opm (not from squid (!))
pos = mean(stim_dpN20_opm.dip.pos,1);
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off

view(12, -10)

%% dmu (analytical)

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\pos275.mat
load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\mom275.mat

%%% mean_residue
Vdata   = mean(stim_dpN20_squid.Vdata,2);
Vmodel  = mean(stim_dpN20_squid.Vmodel,2);
residue = abs(Vdata - Vmodel);

mean_residue_squid = mean(residue); % mean residue across all sensors


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
norm(orient); % check passed

lf = leadfield_squid.leadfield{1,1} * orient;
norm_lf_squid = vecnorm(lf)


%%% calculate dmu
dmu_squid = mean_residue_squid/norm_lf_squid % QUESTION: Is A* cm the units for dmu?
