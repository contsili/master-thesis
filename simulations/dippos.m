% first add Fieltrip to your path: using ft_defaults

%% 1. opm_tangential_radial

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); % headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3) doesn't work since 'cocentric spheres' is not supported as a head model.
sensors = ni2_sensors('type', 'opm_tangential_radial');

% These can be visualized using a few functions from the FieldTrip-toolbox:
figure; hold on;
ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facecolor', 'skin'); alpha 0.5 
ft_plot_sens(sensors, 'elecshape', 'disc', 'elecsize', 1, 'label', 'off');
view([1 1 1])
axis on; grid on
xlabel('x'); ylabel('y'); zlabel('z'); 
 

dippar1 = [0 0 6 1 0 0];

% % Without random sensor noise
% leadfield1 = ni2_leadfield(sensors, headmodel, dippar1);
% ni2_topoplot(sensors, leadfield1); colorbar 

% With random sensor noise
leadfield1 = ni2_leadfield(sensors, headmodel, dippar1, 0.1 * 10^-10);
ni2_topoplot(sensors, leadfield1); colorbar

% a reason why the topoplot
% is not good for opm_tangential_radial is because I have 2 signals from 1
% location. But the lf thinks these are two different locations, so in some
% close locations the signals add up and some others they subtract leading
% to this highly variable topography. The solution would be to add up the
% signals that are picked up by the sensors in the same location:
% leadfield= leadfield2 + leadfield1 and then plot 261 locations, not 522.
% To do that:
% leadfield= leadfield2 + leadfield1;
% ni2_topoplot(sensors1, leadfield); colorbar % opm_tangential + opm_radial
% 
% % Also I can add them up vectorially since the sensors have 90 deg angle
% % difference:
% leadfield= sqrt(leadfield2.^2 + leadfield1.^2);
% ni2_topoplot(sensors1, leadfield); colorbar
% 
% ROBERT SAID THAT THIS IDEA IS FALSE!!!

 
%%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%

data = [];
data.time=linspace(0.01,1,100);
data.avg = repmat(leadfield1,1,100); % keep the same leadfield in 100 time points so I can add a fake time dimension (later I care only about 1 time dimension: cfg.latency = 0.50;). Otherwise I could use ni2_activation. If I do not do that I get the ERROR that the data doesnt represent real "timelock" activity.  
data.label = sensors.label;
data.grad = sensors; % Note: use data.elec for eeg and data.grad for meg
data.dimord = 'chan_time';

cfg = [];
cfg.gridsearch = 'yes';
% cfg.model = 'regional'
cfg.latency = 0.50;
cfg.headmodel = headmodel;
cfg.nonlinear = 'yes';
cfg.numdipoles = 1;

dip = ft_dipolefitting(cfg, data); %this does all the work


ni2_topoplot(sensors,dip.Vdata(:,1)); title('Data');
ni2_topoplot(sensors,dip.Vmodel(:,1)); title('Model');

figure; hold on;
ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facecolor', 'skin'); alpha 0.5 
ft_plot_sens(sensors, 'elecshape', 'disc', 'elecsize', 1, 'label', 'off');
ft_plot_dipole(dip.dip.pos(1,:), mean(dip.dip.mom(1:3,:),2), 'color', 'r')
view([1 1 1])
axis on; grid on
xlabel('x'); ylabel('y'); zlabel('z');  

%%%% Conclusion: the dipole looks at the positive x (like my modeled
%%%% dipole: dippar1 = [0 0 6 1 0 0]), so the fitting works
%%%% even if the topography looked strange, without a clear bipolar pattern

%% 2. opm_tangential

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors1 = ni2_sensors('type', 'opm_tangential');

dippar1 = [0 0 6 1 0 0];

% With random sensor noise
leadfield1 = ni2_leadfield(sensors1, headmodel, dippar1, 0.1 * 10^-10);
ni2_topoplot(sensors1, leadfield1); colorbar % opm_tangential


%%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%

data = [];
data.time=linspace(0.01,1,100);
data.avg = repmat(leadfield1,1,100);
data.label = sensors1.label;
data.grad = sensors1; % Note: use data.elec for eeg and data.grad for meg
data.dimord = 'chan_time';

cfg = [];
cfg.gridsearch = 'yes';
% cfg.model = 'regional'
cfg.latency = 0.50;
cfg.headmodel = headmodel;
cfg.nonlinear = 'yes';
cfg.numdipoles = 1;

dip = ft_dipolefitting(cfg, data); %this does all the work


ni2_topoplot(sensors1,dip.Vdata(:,1)); title('Data');
ni2_topoplot(sensors1,dip.Vmodel(:,1)); title('Model');

figure; hold on;
ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facecolor', 'skin'); alpha 0.5 
ft_plot_sens(sensors1, 'elecshape', 'disc', 'elecsize', 1, 'label', 'off');
ft_plot_dipole(dip.dip.pos(1,:), mean(dip.dip.mom(1:3,:),2), 'color', 'r')
view([1 1 1])
axis on; grid on
xlabel('x'); ylabel('y'); zlabel('z'); 

%% 3. opm_radial

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors2 = ni2_sensors('type', 'opm_radial');

dippar1 = [0 0 6 1 0 0];

% With random sensor noise
leadfield2 = ni2_leadfield(sensors2, headmodel, dippar1, 0.1 * 10^-10);
ni2_topoplot(sensors2, leadfield2); colorbar 


%%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%

data = [];
data.time=linspace(0.01,1,100);
data.avg = repmat(leadfield2,1,100);
data.label = sensors2.label;
data.grad = sensors2; % Note: use data.elec for eeg and data.grad for meg
data.dimord = 'chan_time';

cfg = [];
cfg.gridsearch = 'yes';
% cfg.model = 'regional'
cfg.latency = 0.50;
cfg.headmodel = headmodel;
cfg.nonlinear = 'yes';
cfg.numdipoles = 1;

dip = ft_dipolefitting(cfg, data); %this does all the work


ni2_topoplot(sensors2,dip.Vdata(:,1)); title('Data');
ni2_topoplot(sensors2,dip.Vmodel(:,1)); title('Model');

figure; hold on;
ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facecolor', 'skin'); alpha 0.5 
ft_plot_sens(sensors2, 'elecshape', 'disc', 'elecsize', 1, 'label', 'off');
ft_plot_dipole(dip.dip.pos(1,:), mean(dip.dip.mom(1:3,:),2), 'color', 'r')
view([1 1 1])
axis on; grid on
xlabel('x'); ylabel('y'); zlabel('z'); 


%% 4. ctf275

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors = ni2_sensors('type', 'ctf275');

% These can be visualized using a few functions from the FieldTrip-toolbox:
figure; hold on;
ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facecolor', 'skin'); alpha 0.5 
ft_plot_sens(sensors, 'elecshape', 'disc', 'elecsize', 1, 'label', 'off');
view([1 1 1])
axis on; grid on
xlabel('x'); ylabel('y'); zlabel('z'); 
 
dippar1 = [0 0 6 1 0 0];

% Without random sensor noise
leadfield1 = ni2_leadfield(sensors, headmodel, dippar1);
ni2_topoplot(sensors, leadfield1); colorbar % a reason why the topoplot is not good for opm_tangential_radial is because I have 2 signals from 1 location. But the lf thinks these are two different locations, so in some close locations the signals add up and some others they subtract leading to this highly variable topography. The solution would be to add up the signals that are picked up by the sensors in the same location: leadfield= leadfield2 + leadfield1 and then plot 261 locations, not 522. I do that in the next lines

% % With random sensor noise
% leadfield1 = ni2_leadfield(sensors, headmodel, dippar1, 0.1 * 10^-10);
% ni2_topoplot(sensors, leadfield1); colorbar 



%% TODO: 
% 1. 100 times loop around different nr of channels with Monte Carlo simulations (randomly take out a
% channel IN EVERY ONE of the 20 dipole fits) or mesh_sphere. 
% 2. 100 times loop around different sensor noise
% 3. use qsubselfun
% 4a. For opm do many dipole fits, find σ_opm 
% 4b. For ctf do 1 dipole fit (what sensor noise shall I use?)
% 5. Make the plot that Robert asked (sensors noise - sensors number plot
% and color code with σ_opm/σ_ctf>1 => blue etc...). Note a t-test would
% compare each of the 100*100 σ_opm vs 1 σ_ctf separately and it will give ONLY 1 RESULT (which is if there is
% significant difference or not). We dont want only 1 result!! 

% Maybe I can use Github projects for that.


% Stretch goals: 
% 1. find the typical sensor noise in CTF, 
% 2. add correlatted noise (maybe this is added by adding dipoles inside the head (background brain activity) and outside the head (correlated noise from environmental sources like a car or elevator*))
% * a car or elevator moves regarding the MEG, so it is a moving source (we will model this later)


%% sensors noise - sensors number plot

%% 1 dipole fit for SQUID
sensor_noise_SQUID = 1 * 10^-10; % I picked it randomly
for k = 1:20            
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'ctf275');  

    dippar1 = [0 0 6 1 0 0];
    
    leadfield1 = ni2_leadfield(sensors, headmodel, dippar1, sensor_noise_SQUID);

    %%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%
    data = [];
    data.time=linspace(0.01,1,100);
    data.avg = repmat(leadfield1,1,100); % keep the same leadfield in 100 time points so I can add a fake time dimension (later I care only about 1 time dimension: cfg.latency = 0.50;). Otherwise I could use ni2_activation. If I do not do that I get the ERROR that the data doesnt represent real "timelock" activity.  
    data.label = sensors.label;
    data.grad = sensors; % Note: use data.elec for eeg and data.grad for meg
    data.dimord = 'chan_time';
    
    cfg = [];
    cfg.gridsearch = 'yes';
    % cfg.model = 'regional'
    cfg.latency = 0.50;
    cfg.headmodel = headmodel;
    cfg.nonlinear = 'yes';
    cfg.numdipoles = 1;
    
    dip = ft_dipolefitting(cfg, data); 
    
    dip_location_SQUID = zeros(20,3); % initialise vector
    dip_location_SQUID(k,:) = dip.dip.pos;
end 

sigma_location_SQUID = sqrt((std(dip_location_SQUID(:,1))^2 + std(dip_location_SQUID(:,2))^2 + std(dip_location_SQUID(:,3))^2)  / 3 ); % this formula is taken from Vrba (2000). "Multichannel SQUID biomagnetic systems"
    

%% Many dipole fits for OPM
sensor_noise_OPM = [5 * 10^-10, 7 * 10^-10, 9 * 10^-10, 10 * 10^-10, 11 * 10^-10, 13 * 10^-10, 16 * 10^-10,  20 * 10^-10, 25 * 10^-10];
sensor_number_OPM = 21:15:201;
sigma_location_OPM = zeros(length(sensor_noise_OPM), length(sensor_number_OPM)); % initialise vector
sigma_location_relative = zeros(length(sensor_noise_OPM), length(sensor_number_OPM)); % initialise vector

for i = 1:length(sensor_noise_OPM)
    for j = 1:length(sensor_number_OPM)      
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            fprintf('\n\n')
            fprintf('Sensor Number (OPM): %d\n', sensor_number_OPM(j));
            fprintf('Sensor Noise (OPM): %d\n', i);
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            fprintf('\n\n')

            headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
            sensors = ni2_sensors('type', 'opm_tangential_radial', 'sensor_number', sensor_number_OPM(j));
             
            dippar1 = [0 0 6 1 0 0];
            
            leadfield1 = ni2_leadfield(sensors, headmodel, dippar1);

            %%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%
            data = [];
            data.time=linspace(0.01,1,100); % 100 dipole fit iterations
            data.avg = repmat(leadfield1,1,length(data.time))+ sensor_noise_OPM(i) * randn(1, length(data.time)); % keep the same leadfield in 100 time points so I can add a fake time dimension (later I care only about 1 time dimension: cfg.latency = 0.50;). Otherwise I could use ni2_activation. If I do not do that I get the ERROR that the data doesnt represent real "timelock" activity.  
            data.label = sensors.label;
            data.grad = sensors; % Note: use data.elec for eeg and data.grad for meg
            data.dimord = 'chan_time';
            
            cfg = [];
            cfg.gridsearch = 'no';
            cfg.dip.pos = dippar1(1:3);
            cfg.model = 'moving';
            cfg.latency = 'all';
            cfg.headmodel = headmodel;
            cfg.nonlinear = 'yes';
            cfg.numdipoles = 1;
%             cfg.feedback = 'no'; % Purpose: do not print the results of
%             ft_dipolefitting in the command window. Result: It does not
%             work
            cfg.showcallinfo = 'no'; % the time and memory is not printed
            
            dip = ft_dipolefitting(cfg, data); 

            dip_location_OPM = zeros(length(data.time),3); % initialise vector
            from_structure_to_vector = extractfield(dip.dip, 'pos');
            dip_location_OPM = reshape(from_structure_to_vector, [3,length(data.time)])';
         
        
            sigma_location_OPM(i,j) = sqrt((std(dip_location_OPM(:,1))^2 + std(dip_location_OPM(:,2))^2 + std(dip_location_OPM(:,3))^2)  / 3 ); % this formula is taken from Vrba (2000). "Multichannel SQUID biomagnetic systems"
            sigma_location_relative(i,j) = sigma_location_OPM(i,j) / sigma_location_SQUID;
            % sigma_location_log = log(sigma_location_OPM / sigma_location_SQUID);
    
            if sigma_location_relative(i,j) > 2 % we care about the "green" line, so when sigma_location_relative(i,j) = 1, we do not care about too high values. I chose sigma_location_relative(i,j) > 2 to get an equidistant range around 1 (from 0 to 1 OPM wins, from 1 to 2 SQUID wins)
                sigma_location_relative(i,j) = nan;
            end

    end 
end



%% plot 

% Plot1: sensor noise OPM - sensor number OPM - sigma
subplot(2,length(sensor_noise_OPM),1:length(sensor_noise_OPM))
h = imagesc(sensor_number_OPM, sensor_noise_OPM, sigma_location_relative);
colormap('jet'); 
set(gca,'YDir','normal') % flip the y-axis to be in ascending order

% Set the axis labels
xlabel('sensor number OPM');
ylabel('sensor noise OPM');
title('Dipole fit uncertainty for OPM vs SQUID')

% Add a colorbar with the label "sigma_OPM/sigma_SQUID"
c = colorbar;
c.Label.Interpreter = 'tex';
c.Label.String = '\sigma_{OPM}/\sigma_{SQUID}';

% Modify the colorbar to show "SQUID wins" when the color is blue and "OPM wins" when the color is red
colormap(c, 'jet');
c.Ticks = [0.6, 1.75];
c.TickLabels = {'OPM wins', 'SQUID wins'};

% make NaNs totally transparent (i.e., make them white)
set(h, 'AlphaData', ~isnan(sigma_location_relative));



% Plot2: sigma - sensor number OPM
for i=1:length(sensor_noise_OPM)
    subplot(2,length(sensor_noise_OPM),length(sensor_noise_OPM)+i)
    plot(sensor_number_OPM, sigma_location_relative(i,:),'b*')
    ylabel('\sigma_{OPM}/\sigma_{SQUID}')
    xlabel('sensor number OPM')

    h = sprintf('Sensor noise: %g', sensor_noise_OPM(i)); % note that fprintf() does not work.
    title(h);

    % Fit a linear model
    y = sigma_location_relative(i, :);
    x = sensor_number_OPM;
    
    % Fit a linear model (y = mx + b)
    p = polyfit(x(~isnan(y)), y(~isnan(y)), 4);
    
    % Calculate the fitted values
    fittedValues = polyval(p, x);
    
    % Plot the fitted line
    hold on;
    plot(x, fittedValues, 'r', 'LineWidth', 2);
    

%     % Add the fitted equation as text below the plot
%     fittedEquation = sprintf('y = %.4fx + %.4f', p(1), p(2));
% 
%     % Get current axes
%     ax = gca; 
%     
%     % Add text on the right upper corner of each subplot
%     text(ax.XLim(2), ax.YLim(2), fittedEquation, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
end


