function sens = ni2_sensors(varargin)

% NI2_SENSORS create a fictive sensor-array.
%
% Use as
%  sens = ni2_sensors('type', senstype);
%
% Where senstype can be any of the following:
%  'eeg' generates a 91-channel eeg sensor array
% 
%  'meg' generates a 301-channel meg magnetometer sensor array
% 
%  'opm_radial' generates a 301-channel opm-meg magnetometer sensor array with
%  one radial channel for each sensor (the difference from 'meg' is that
%  'opm_radial' are closer to the scalp)
% 
%  'opm_tangential' generates a 301-channel opm-meg magnetometer sensor array with
%  one tangential channel for each sensor
% 
% 'opm_tangential_radial' generates a 301-channel opm-meg magnetometer sensor array with
%  two channels for each sensor, one tangential and one radial (like the
%  FieldLine system we bought in DCCN)


type   = ft_getopt(varargin, 'type', 'eeg');
jitter = ft_getopt(varargin, 'jitter', 0);
n      = ft_getopt(varargin, 'n', 162); % number of vertices for the sphere, determines the number of electrodes, this is the old default
NumberOfSensors = ft_getopt(varargin, 'sensor_number', 0);

switch type
  case 'eeg'
    % create an electrode array
    [chanpos, tri] = mesh_sphere(n);
    chanpos        = chanpos*10;
    chanpos(chanpos(:,3)<0,:) = [];

    if jitter
      [th,phi,r] = cart2sph(chanpos(:,1), chanpos(:,2), chanpos(:,3));
      shift1 = 2*jitter*(rand(numel(th),1)  - 1);
      shift2 = 2*jitter*(rand(numel(phi),1) - 1);
      [chanpos(:,1),chanpos(:,2),chanpos(:,3)] = sph2cart(th+shift1(:),phi+shift2(:),r);
    end

    sens.unit    = 'cm';
    sens.coordsys = 'neuromag';
    sens.chanpos = chanpos;
    sens.elecpos = chanpos;
    for k = 1:size(chanpos,1)
      sens.label{k,1}    = sprintf('eeg% 02d', k);
      sens.chantype{k,1} = 'eeg';
    end
    sens = ft_datatype_sens(sens);
    

    case {'meg'}
    % create a magnetometer array
    [chanpos, tri] = icosahedron642;
     coilori        = chanpos;
    
    % rotated_points now contains the rotated points
    chanpos        = chanpos*12.5; % why *12.5? Is there a mathematical formula for that?
    chanpos(:,3)   = chanpos(:,3) - 1.5;

    z = chanpos(:,3);
    coilori        = coilori(z>0,:);
    chanpos        = chanpos(z>0,:);
    nchan          = size(chanpos,1);

    sens.unit    = 'cm';
    sens.coordsys = 'neuromag';
    sens.chanpos = chanpos;
    sens.chanori = coilori;
    sens.coilpos = chanpos;
    sens.coilori = coilori;
    sens.tra     = eye(nchan);
    for k = 1:nchan
      sens.label{k,1}    = sprintf('meg% 03d', k);
      sens.chantype{k,1} = 'megmag';
    end
    sens = ft_datatype_sens(sens);


    case {'opm_tangential_radial'}
    % create a magnetometer array
    [chanpos, tri] = icosahedron642;
    
    % Keep only the upper semi-sphere
    z = chanpos(:,3);
    chanpos       = chanpos(z>0,:);
    
    % Manually determine the number of sensors (if NumberOfSensors = 32
    % then number of channels will be 64, which means chanpos: 64 x 3, and not 32 x 3)
    if NumberOfSensors
        % Randomly take out a sensors (Monte-Carlo simulation)
        randomIndices = randperm(size(chanpos, 1), NumberOfSensors);
        chanpos       = chanpos(randomIndices,:);
    end
    
    coilori = zeros(2*size(chanpos,1),3);
    coilori(1:size(chanpos,1),:)        = chanpos;
    angle =pi/2;
    j=1;
    for i = size(chanpos, 1)+1 : 2*size(chanpos, 1) 
        % Convert to spherical coordinates
        [phi, theta, rho] = cart2sph(chanpos(j, 1), chanpos(j, 2), chanpos(j, 3));
        % Define the rotation axis as the phi vector
        rotation_axis = [sin(phi) -cos(phi) 0];
        % Define the rotation matrix using Rodrigues' rotation formula
        R = eye(3) + sin(angle) * cross_matrix(rotation_axis) + (1 - cos(angle)) * cross_matrix(rotation_axis)^2;
        % Rotate the point
        coilori(i, :) = (R * chanpos(j, :)')';
        j=j+1;
    end
    
    % rotated_points now contains the rotated points
    chanpos        = chanpos*10.5; % I choose *10 as in the eeg, because opm are close to the head

    chanpos=[chanpos;chanpos];    
    nchan          = size(chanpos,1);    

    sens.unit    = 'cm';
    sens.coordsys = 'neuromag';
    sens.chanpos = chanpos;
    sens.chanori = coilori;
    sens.coilpos = chanpos;
    sens.coilori = coilori;
    sens.tra     = eye(nchan);
    for k = 1:nchan
      sens.label{k,1}    = sprintf('meg% 03d', k);
      sens.chantype{k,1} = 'megmag';
    end
    sens = ft_datatype_sens(sens);


    case {'opm_tangential'}
    % create a magnetometer array
    [chanpos, tri] = icosahedron642;

    coilori = zeros(size(chanpos));
    angle =pi/2;
    for i = 1:size(chanpos, 1)
        % Convert to spherical coordinates
        [phi, theta, rho] = cart2sph(chanpos(i, 1), chanpos(i, 2), chanpos(i, 3));
        % Define the rotation axis as the phi vector
        rotation_axis = [sin(phi) -cos(phi) 0];
        % Define the rotation matrix using Rodrigues' rotation formula
        R = eye(3) + sin(angle) * cross_matrix(rotation_axis) + (1 - cos(angle)) * cross_matrix(rotation_axis)^2;
        % Rotate the point
        coilori(i, :) = (R * chanpos(i, :)')';
    end
    
    % rotated_points now contains the rotated points
    chanpos        = chanpos*10.5; % I choose *10 as in the eeg, because opm are close to the head

    z = chanpos(:,3);
    coilori        = coilori(z>0,:);
    chanpos        = chanpos(z>0,:);
    nchan          = size(chanpos,1);

    sens.unit    = 'cm';
    sens.coordsys = 'neuromag';
    sens.chanpos = chanpos;
    sens.chanori = coilori;
    sens.coilpos = chanpos;
    sens.coilori = coilori;
    sens.tra     = eye(nchan);
    for k = 1:nchan
      sens.label{k,1}    = sprintf('meg% 03d', k);
      sens.chantype{k,1} = 'megmag';
    end
    sens = ft_datatype_sens(sens);


     case {'opm_radial'}
    % create a magnetometer array
    [chanpos, tri] = icosahedron642;
     coilori        = chanpos;
    
    % rotated_points now contains the rotated points
    chanpos        = chanpos*10;

    z = chanpos(:,3);
    coilori        = coilori(z>0,:);
    chanpos        = chanpos(z>0,:);
    nchan          = size(chanpos,1);

    sens.unit    = 'cm';
    sens.coordsys = 'neuromag';
    sens.chanpos = chanpos;
    sens.chanori = coilori;
    sens.coilpos = chanpos;
    sens.coilori = coilori;
    sens.tra     = eye(nchan);
    for k = 1:nchan
      sens.label{k,1}    = sprintf('meg% 03d', k);
      sens.chantype{k,1} = 'megmag';
    end
    sens = ft_datatype_sens(sens);


  case 'ctf151'
    load('ctf151');
    % only keep the normal channels

    sel = strcmp(sens.chantype, 'meggrad');
    sens.label    = sens.label(sel);
    sens.chantype = sens.chantype(sel);
    sens.chanunit = sens.chanunit(sel);
    sens.tra      = sens.tra(sel,:);
    sens.chanpos  = sens.chanpos(sel,:);
    sens.chanori  = sens.chanori(sel,:);


  case 'ctf275'
    load('ctf275');

    % only keep the normal MEG channels
    sel = strcmp(sens.chantype, 'meggrad');
    sens.label    = sens.label(sel);
    sens.chantype = sens.chantype(sel);
    sens.chanunit = sens.chanunit(sel);
    sens.tra      = sens.tra(sel,:);
    sens.chanpos  = sens.chanpos(sel,:);
    sens.chanori  = sens.chanori(sel,:);



  case 'bti248'
    load('bti248');

    % only keep the normal MEG channels
    sel = strcmp(sens.chantype, 'megmag');
    sens.label    = sens.label(sel);
    sens.chantype = sens.chantype(sel);
    sens.chanunit = sens.chanunit(sel);
    sens.tra      = sens.tra(sel,:);
    sens.chanpos  = sens.chanpos(sel,:);
    sens.chanori  = sens.chanori(sel,:);

  case 'neuromag306'
    load('neuromag306');

  case 'eeg1020'
    load('eeg1020');

  case 'eeg1010'
    load('eeg1010');

  otherwise
    error('unsupported type % s', type);
end


    function C = cross_matrix(v) % used in the Rodrigues formula
    % This function returns the cross product matrix of a vector v
    C = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end
end
