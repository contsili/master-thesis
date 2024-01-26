beta2 = ft_read_json('beta2.geo');

% Here is the beta2 helmet nominal values. The 12 parameters are in the MEGIN
% convention of specifying position in m (3 parameters) followed by the directions of
% the 3 orthogonal axis for orientation (9 parameters, ex (unit vector), ey (unit
% vector), ez (unit vector))
%
% Note that the we dropped the “dash” in the “hole names” in this file. The channel
% names have the following convention: <hole name>_<axis>-<chassis location>
%
% Thus, in the beta2.geo file, L1-01 becomes L111 and if you are measuring the bz
% direction, -> L111_bz. If L111 is in plugged into the first port and first card in
% the chassis, the full channel name is L111_bz-s1


%%

label = fieldnames(beta2.location);

clear coil*

for i=1:numel(label)

  coilpos(i,:) = beta2.location.(label{i})(2:4);
  coilorix(i,:) = beta2.location.(label{i})(5:7);
  coiloriy(i,:) = beta2.location.(label{i})(8:10);
  coiloriz(i,:) = beta2.location.(label{i})(11:13);
end

%%

gradx = [];
gradx.coilpos = coilpos;
gradx.coilori = coilorix;
gradx.tra     = eye(144);
gradx.label = {};
for i=1:144
  gradx.label{i} = [label{i}([1 2 4 5]) '_bx'];
end

grady = [];
grady.coilpos = coilpos;
grady.coilori = coiloriy;
grady.tra     = eye(144);
grady.label = {};
for i=1:144
  grady.label{i} = [label{i}([1 2 4 5]) '_by'];
end

gradz = [];
gradz.coilpos = coilpos;
gradz.coilori = coiloriz;
gradz.tra     = eye(144);
gradz.label   = {};
for i=1:144
  gradz.label{i} = [label{i}([1 2 4 5]) '_bz'];
end

for i=1:144
  gradx.coilori(i,:) = gradx.coilori(i,:) / norm(gradx.coilori(i,:));
  grady.coilori(i,:) = grady.coilori(i,:) / norm(grady.coilori(i,:));
  gradz.coilori(i,:) = gradz.coilori(i,:) / norm(gradz.coilori(i,:));
end

%%

figure; ft_plot_sens(gradx, 'orientation', 1, 'label', 'label');
figure; ft_plot_sens(grady, 'orientation', 1, 'label', 'label');
figure; ft_plot_sens(gradz, 'orientation', 1, 'label', 'label');

%%

cfg = [];
grad = ft_appendsens(cfg, gradx, grady, gradz);

%%

figure; ft_plot_sens(grad, 'orientation', 1, 'label', 'label');

%%

save gradx gradx
save grady grady
save gradz gradz
save grad  grad

%%


cfg = [];
cfg.grad = gradx;
cfg.skipcomnt = 'yes';
cfg.skipscale = 'yes';
cfg.outline = 'helmet';
cfg.mask = 'helmet';

layoutx = ft_prepare_layout(cfg);
layoutx.pos(:,2) = 1.2*layoutx.pos(:,2)+0.3;

layouty = ft_prepare_layout(cfg);
layouty.pos(:,2) = 1.2*layouty.pos(:,2)+0.3;

layoutz = ft_prepare_layout(cfg);
layoutz.pos(:,2) = 1.2*layoutz.pos(:,2)+0.3;

figure; ft_plot_layout(layoutz)

%%

cfg = [];
cfg.distance = 1;
layout = ft_appendlayout(cfg, layoutx, layouty, layoutz);

figure; ft_plot_layout(layout)

%%

save layoutx layoutx
save layouty layouty
save layoutz layoutz
save layout  layout

