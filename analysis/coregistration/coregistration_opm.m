%% CO-REGISTRATION

addpath('H:/common/matlab/fieldtrip')

ft_defaults

%% Read the 3D scan, assign a meaningful coordinate system, and erase the irrelevant parts

% Read the 3D scan

% 001-2 scan
scan      = ft_read_headshape('M:\Documents\recordings\opm-MN\optical-scan\001-02\Model\Model.obj');
scan.unit = 'm'; % the estimated 'dm' is not correct

figure; hold on;
ft_plot_headshape(scan);
ft_plot_axes(scan);
lighting gouraud
material dull
light

% % 001 scan
% scan      = ft_read_headshape('M:\Documents\recordings\opm-MN\optical-scan\001\Model\Model.obj');
% scan.unit = 'm'; % the estimated 'dm' is not correct
% 
% figure; hold on;
% ft_plot_headshape(scan);
% ft_plot_axes(scan);
% lighting gouraud
% material dull
% light

% CONCLUSION: I will keep scan 001-2 since: 
% 1. it is taken in front of the subject face while scan 001 is taken from the side. That way, I can localise LPA better. 
% 2. The subject has a neutral face reaction and
% the rim of the helmet is clearer.


% assign a meaningful coordinate system: approximately align the mesh to a RAS coordinate system,
% by clicking on 'dummy' nas/lpa/rpa
cfg          = [];
cfg.method   = 'fiducial';
cfg.coordsys = 'neuromag'; % RAS (Right, Anterior, Superior) and hence we specify cfg.coordsys='neuromag'
scan_aligned = ft_meshrealign(cfg, scan);

figure; hold on;
ft_plot_headshape(scan_aligned);
ft_plot_axes(scan_aligned);
view([125 10]);
lighting gouraud
material dull
light


% erase the irrelevant parts (erase body & cables)
cfg         = [];
cfg.method  = 'box';
cfg.selection = 'inside';
scan_head   = ft_defacemesh(cfg, scan_aligned); 

figure; hold on;
ft_plot_headshape(scan_head);
ft_plot_axes(scan_head);
view([125 10]);
lighting gouraud
material dull
light



% get the face 
cfg        = [];
cfg.method = 'box';
cfg.selection = 'inside';
scan_face  = ft_defacemesh(cfg, scan_head); 

figure; hold on;
ft_plot_headshape(scan_face);
ft_plot_axes(scan_face);
view([125 10]);
lighting gouraud
material dull
light


% get the helmet
cfg           = [];
cfg.method    = 'box';
cfg.selection = 'outside';
scan_helmet   = ft_defacemesh(cfg, scan_head);

figure; hold on;
ft_plot_headshape(scan_helmet);
ft_plot_axes(scan_helmet);
view([125 10]);
lighting gouraud
material dull
light

% helmet still has some face parts

%% Read the anatomical MRI, assign a well defined head coordinate system

% 1. Read the anatomical MRI
mri          = ft_read_mri('M:\Documents\recordings\squid-MN\dicoms\00001_1.3.12.2.1107.5.2.19.45416.2022110716263882460203497.IMA');

% check the coordinate system. I expect that the origin of the coordinate system is ill-defined,
% as that depends how the subject was lying in the MRI scanner and how the scanned volume was configured
ft_determine_coordsys(mri, 'interactive', 'no'); 
view([0 100]); % It is RAS (Right, Anterior, Superior). The NAS is not defined well 


% 2. assign a well defined head coordinate system
% Now I define the coordinate system manually by clicking at NAS, RPA, LPA.
cfg = [];
cfg.coordsys = 'neuromag'; % in the squid analysis I used 'ctf'!
mri_realigned = ft_volumerealign(cfg, mri);

% check the coordinate system
ft_determine_coordsys(mri_realigned, 'interactive', 'no'); % Now the coordinate system origin is more inferiorly than before

%% Interactive alignment of the face - extracted from the 3D scan - with the face extracted from the anatomical MRI

% segment the scalp
cfg          = [];
cfg.output   = 'scalp';
seg          = ft_volumesegment(cfg, mri_realigned);

% create a mesh for the scalp
cfg             = [];
cfg.tissue      = 'scalp';
cfg.numvertices = 10000;
mri_face        = ft_prepare_mesh(cfg, seg);
mri_face        = ft_convert_units(mri_face, 'm');


% align the 3D scan face to the face extracted from the anatomical mri
cfg             = [];
cfg.method      = 'interactive';
cfg.headshape   = mri_face;
cfg.meshstyle   = {'edgecolor', 'k', 'facecolor', 'skin'};
scan_face_aligned = ft_meshrealign(cfg, scan_face); % do not look so much on the whole nose. Look to fit well the nasion part

figure; hold on;
ft_plot_headshape(mri_face, 'facealpha', 0.4);
ft_plot_mesh(scan_face_aligned, 'facecolor','skin');
view([125 10]);
lighting gouraud
material dull
light

%% Interactive alignment of the helmet - extracted from the 3D scan - with the reference helmet and sensors (rim of the helmet)

helmet_rim = ft_read_headshape('fieldlinebeta2_helmet_rim.mat');
helmet_rim.coordsys = 'ras';

cfg = [];
cfg.method = 'interactive';
cfg.headshape = helmet_rim;
cfg.meshstyle = {'edgecolor', 'none', 'facecolor', [1 0.5 0.5]};
scan_helmet_aligned = ft_meshrealign(cfg, scan_helmet);

figure; hold on;
ft_plot_mesh(helmet_rim, 'edgecolor', 'none', 'facecolor', [0.5 0.5 1], 'facealpha', 0.4);
ft_plot_mesh(scan_helmet_aligned, 'edgecolor', 'none', 'facecolor', [1 0.5 0.5], 'facealpha', 0.4);
view([145 10]);
lighting gouraud
material dull
light

%% calculate and apply the transformation that aligns the helmet_rim with the anatomical MRI.

% calculate the transformation
transform_scan2helmet = scan_helmet_aligned.cfg.transform;
transform_scan2face   = scan_face_aligned.cfg.transform;
transform_helmet2face = transform_scan2face/transform_scan2helmet; 
% EXPLAIN THE MATH:
% 1/transform_scan2helmet: helmet -> scan, 
% transform_scan2face: scan -> face
% transform_scan2face/transform_scan2helmet: helmet -> scan -> face


% apply the transformation
fieldlinebeta2 = ft_read_sens('20240305_134446_sub-001_file-MedianNerve02_raw.fif');
fieldlinebeta2.coordsys = 'ras';

fieldlinebeta2_head = ft_transform_geometry(transform_helmet2face, fieldlinebeta2);

figure; hold on;
ft_plot_sens(fieldlinebeta2_head);
ft_plot_headshape(mri_face, 'facecolor', [0.5 0.5 1], 'facealpha', 0.4, 'edgecolor', 'none');
view([125 10]);
lighting gouraud
material dull
light

% PROBLEM: run1: left temporal and parieto-occipital sensors appear to be
% in the scalp


%% append sensors

% run1
sens1 = ft_read_sens('20240305_133110_sub-001_file-MedianNerve01_raw.fif');

% run2
sens2 = ft_read_sens('20240305_134446_sub-001_file-MedianNerve02_raw.fif');
% sens2.label(find(strcmp(sens2.label, {'R105_bz'}))) = []; % delete 'R105_bz' which is a duplicate channel
% 
% PROBLEM: i need to delete 'R105_bz' also from sens3.chanori etc... A
% solution could be (but this function does not exist):
% cfg =[]; 
% cfg.channel ={'all', '-R105_bz', '-R603_bz'} 
% sens2 = ft_selectsens(cfg. sens2)
%
% For faster I do it at the ft_dipolefitting()

% run3
sens3 = ft_read_sens('20240305_135949_sub-001_file-MedianNerve03_raw.fif');

% run4
sens4 = ft_read_sens('20240305_141224_sub-001_file-MedianNerve004_raw.fif');
% sens4.label(find(strcmp(sens4.label, {'R603_bz'}))) = []; % delete 'R603_bz' which is a duplicate channel

% run5
sens5 = ft_read_sens('20240305_142537_sub-001_file-MedianNerve05_raw.fif');

% run6
sens6          = ft_read_sens('20240305_143835_sub-001_file-MedianNerve06_raw.fif');
sens6.label(6) = {'L101_bz'};

% append sensors. This function omits automatically one of the 2 duplicate channels.
cfg      = [];
combined = ft_appendsens(cfg, sens1, sens2, sens3, sens4, sens5, sens6);
combined.coordsys = 'ras';

% CONCLUSION: 5 out of the 144 sensors missing as expected. '-R105_bz',
% '-R603_bz' will be taken care of at the ft_dipolefitting. So in total 7 out of the 144 sensors missing

%% plot the whole head sensor coverage

fieldlinebeta2_head = ft_transform_geometry(transform_helmet2face, combined);
fieldlinebeta2_head = ft_convert_units(fieldlinebeta2_head, 'cm');

figure; hold on;
ft_plot_sens(fieldlinebeta2_head, 'label', 'on');
ft_plot_headshape(mri_face, 'facecolor', [0.5 0.5 1], 'facealpha', 0.4, 'edgecolor', 'none');
view([125 10]);
lighting gouraud
material dull
light

% CONCLUSION: L101 and L214 appear too deep. I am going to ignore them for
% now. They both seem to be problem of the HED software helmetscan

% Note: L214 was probably mixed with L312 by the software at run2. I can
% still use L214 from run3 (but I need to do that before I run ft_appendsens())




%% HEADMODEL 

cfg                = [];
cfg.output         = 'brain';
% cfg.scalpthreshold = 0.25;
mri_segmented          = ft_volumesegment(cfg, mri_realigned);

% copy the anatomy into the segmented mri. That way we can use sourceplot
% after to plot the 3 segmentations on top of the anatomical mri
mri_segmented.anatomy = mri_realigned.anatomy;

mri_segmented = ft_convert_units(mri_segmented, 'cm');

cfg = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, mri_segmented);


cfg                = [];
cfg.tissue         = {'brain'};
mesh_brain          = ft_prepare_mesh(cfg, mri_segmented);


cfg        = [];
cfg.method = 'singlesphere';
headmodel_sphere  = ft_prepare_headmodel(cfg, mesh_brain);

headmodel_sphere = ft_convert_units(headmodel_sphere, 'cm');

figure;
hold on
ft_plot_headmodel(headmodel_sphere, 'facealpha', 0.5)
ft_plot_sens(fieldlinebeta2_head)
camlight
% CONCLUSION: Superior sensors too close to the brain.





%% DIPOLE FIT

append_tl.avg = append_tl.trial;

cfg            = [];
cfg.latency    = [0.035 0.05];
cfg.unit       = 'cm';
cfg.gridsearch = 'yes';
cfg.grad       = fieldlinebeta2_head;
cfg.channel    = {'all', '-R105_bz', '-R603_bz', '-L101_bz', '-L214_bz'}; % Take care of the duplicate channels for run2 and run4 & the channels that appear closer to the brain than expected. In the end I am using 144-9=135 channels
cfg.headmodel  = headmodel_sphere;
stim_dpN20_opm = ft_dipolefitting(cfg, append_tl);

hold on

ft_plot_dipole(stim_dpN20_opm.dip.pos, mean(stim_dpN20_opm.dip.mom(1:3,:),2), 'color', 'b')

pos = mean(stim_dpN20_opm.dip.pos,1);
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_segmented.anatomy, 'transform', mri_segmented.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)

ft_plot_crosshair(pos, 'color', [1 1 1]/2);

axis tight
axis off

view(12, -10)

% CONCLUSION: OPM dipole pointing to the opposite direction is expected as the
% colors in the topo for the OPMs were reversed. OPM dipole appears a bit
% more posteriorly than SQUID

% in topoplotER I used 144-7 sensors while in the dipolefitting I used
% 144-9 sensors (since '-L101_bz' and '-L214_bz' were too deep)

%% find the modelled topo based on the dipole I fitted. Compare it with measured topo

%%% model (Vmodel) 
% way1
cfg =[];
cfg.sourcemodel.pos = pos;
cfg.grad            = fieldlinebeta2_head;
cfg.channel         = {'all', '-R105_bz', '-R603_bz', '-L101_bz', '-L214_bz'}; 
cfg.headmodel       = headmodel_sphere;
leadfield           = ft_prepare_leadfield(cfg);

Vmodel = leadfield.leadfield{1,1}*stim_dpN20_opm.dip.mom;
Vmodel = mean(Vmodel, 2);

cfg         = [];
cfg.layout  = 'fieldlinebeta2bz_helmet';
cfg.channel = {'all', '-R105_bz', '-R603_bz', '-L101_bz', '-L214_bz'}; 

data        = [];
data.avg    = Vmodel;
data.label  = stim_dpN20_opm.label;
data.time   = 1;
data.dimord = 'chan_time';

ft_topoplotER(cfg, data)

% way2
Vmodel = mean(stim_dpN20_opm.Vmodel,2);

cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.channel         = {'all', '-R105_bz', '-R603_bz', '-L101_bz', '-L214_bz'}; 

data        = [];
data.avg    = Vmodel;
data.label  = stim_dpN20_opm.label;
data.time   = 1;
data.dimord = 'chan_time';

ft_topoplotER(cfg, data)



%%% data (measured topo - Vdata)
cfg        = [];
cfg.xlim   = [0.035 0.050];
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.zlim   = 'maxabs';
ft_topoplotER(cfg, append_tl); 


%%% residue
Vdata   = mean(stim_dpN20_opm.Vdata,2);
Vmodel  = mean(stim_dpN20_opm.Vmodel,2);
residue = abs(Vdata - Vmodel);

cfg = [];
cfg.layout = 'fieldlinebeta2bz_helmet';
cfg.channel         = {'all', '-R105_bz', '-R603_bz', '-L101_bz', '-L214_bz'}; 

data        = [];
data.avg    = residue;
data.label  = stim_dpN20_opm.label;
data.time   = 1;
data.dimord = 'chan_time';

ft_topoplotER(cfg, data)


%% dmu (analytical)

%%% mean_residue
Vdata   = mean(stim_dpN20_opm.Vdata,2);
Vmodel  = mean(stim_dpN20_opm.Vmodel,2);
residue = abs(Vdata - Vmodel);

mean_residue = mean(residue); % mean residue across all sensors


%%% norm_lf
cfg =[];
cfg.sourcemodel.pos = pos;
cfg.grad            = fieldlinebeta2_head;
cfg.channel         = {'all', '-R105_bz', '-R603_bz', '-L101_bz', '-L214_bz'}; 
cfg.headmodel       = headmodel_sphere;
leadfield           = ft_prepare_leadfield(cfg);

% QUESTION: leadfield.leadfield{1,1} is a (nr of sensors)x3 matrix. How am I
% going to calculate 1 number for the norm and integrate the x,y and z orientations of the leadfield?

mom = mean(stim_dpN20_opm.dip.mom  , 2);
norm_mom = norm(mom);

orient = mom ./ norm_mom;
norm(orient); % check passed

lf = leadfield.leadfield{1,1} * orient;
norm_lf = vecnorm(lf)


%%% calculate dmu
dmu = mean_residue/norm_lf % QUESTION: Is A* cm the units for dmu?

% dmu =
% 
%    7.6652e-06