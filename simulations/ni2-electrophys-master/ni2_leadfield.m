function leadfield = ni2_leadfield(sens, headmodel, dippar, sigma)
% NI2_LEADFIELD generates a leadfield for a given source (or set of
% sources) using a specified sensor array and volume conductor model.
%
% Use as:
%   leadfield = ni2_leadfield(sens, headmodel, dippar, sigma)
%
% Input arguments:
%   - sens      = a sensor array, obtained with ni2_sensors
%   - headmodel = a volume conductor model, obtained with ni2_headmodel
%   - dippar    = Nx6, or Nx3 matrix with dipole parameters.
%   - sigma     = Standard deviation of sensor noise (optional)
%
% Each row in the dippar-matrix represents a source, the first 3 columns
% are the position parameters (x,y,z coordinates in Cartesian space), and
% the options 4-6 columns are the dipole moment parameters (x,y,z
% orientation and amplitude).
% 
% Output argument:
%   - leadfield = MxN or Mx(Nx3) matrix with the leadfield, where M is the
%                 number of sensors, and N the number of sources

ndip = size(dippar, 1); %number of dipoles

cfg = [];
if strcmp(sens.chantype{1}, 'eeg')
  cfg.elec = sens;
else
  cfg.grad = sens;
end
cfg.headmodel   = headmodel;
cfg.sourcemodel.pos    = dippar(:,1:3);
cfg.sourcemodel.inside = 1:ndip;
cfg.sourcemodel.unit   = headmodel.unit;
cfg.reducerank  = 'no';
if strcmp(headmodel.type, 'singleshell')
  cfg.singleshell.batchsize = 2500;
end
leadf           = ft_prepare_leadfield(cfg);

if size(dippar, 2) == 6
  leadfield = zeros(size(leadf.leadfield{1}, 1), ndip);
  for k = 1:ndip
    leadfield(:, k) = leadf.leadfield{k} * dippar(k, 4:6)';

    % Add Gaussian sensor noise if sigma is provided
    if nargin == 4
      leadfield(:, k) = leadfield(:, k) + sigma * randn(size(leadfield, 1), 1);
    end

  end
else
  leadfield = zeros(size(leadf.leadfield{1}, 1), ndip * 3);
  for k = 1:ndip
    indx = (k-1) * 3 + (1:3);
    leadfield(:, indx) = leadf.leadfield{k};
  end
end
