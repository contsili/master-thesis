%% Analytical solution of sensor number when we assume a certain noise for the OPMs

% This was my try to use eq. (4) (moment uncertainty opm = moment
% uncertainty squid -> both systems have the SAME performance) and find a 
% specific number of opm sensors. 
% 
% In the end:
% 1. this was not possible as ||L|| - number
% of sensors cannot be modeled with an equation. 
% 
% 2. and was not informative since it depends on where we place the dipole. 
% Ashley does not know where the dipole is. The most important is to get a 
% general feeling of the relation between σ, sensor noise and sensor number 
% to get future strategies for the OPMs. 


%% 13/11/2013
% 1. To find the sensor_number OPM I do: sensor_number_OPM = leadfield_pinv_OPM / leadfield_pinv_OPM_default_mean;. So I divide with the mean(leadfield)

%% Way1: Use pinv

sensor_noise_OPM = [1e-10, 5e-10, 7e-10, 9e-10, 1e-9, 1.1e-9, 1.3e-9, 1.6e-9, 2e-9, 2.5e-9]; % Fieldline paper: the remaining average noise floor (after SSP, look figure 5B) is 20 ± 5 fT/Hz^1/2.
sensor_noise_SQUID =  1e-10; % sensor_noise_SQUID =  5 * 10^-15; % In intranet it says: The noise levels of MEG sensors should be below 10 fT, usually around 4-7fT. (https://intranet.donders.ru.nl/index.php?id=meg-lab-procedure-noiserecording&no_cache=1&sword_list%5B%5D=lab)


headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors = ni2_sensors('type', 'ctf275');  
dippar = [0 0 6 1 0 0];
leadfield = ni2_leadfield(sensors, headmodel, dippar, sensor_noise_SQUID);
leadfield_pinv_SQUID = pinv(leadfield);

for i = 1:length(sensor_noise_OPM)
    % Idea from Brookes et al. (2021) which is the same as eq. (4) (see photo in https://github.com/contsili/master-thesis/issues/2#issuecomment-1799608261)
    leadfield_pinv_OPM(i) = norm(leadfield_pinv_SQUID) * sensor_noise_SQUID / sensor_noise_OPM(i);
    
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); % headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3) doesn't work since 'cocentric spheres' is not supported as a head model.
    sensors = ni2_sensors('type', 'opm_tangential_radial');
    dippar1 = [0 0 6 1 0 0]; 
    leadfield1 = ni2_leadfield(sensors, headmodel, dippar1, sensor_noise_OPM(i));
    
    % % These can be visualized using a few functions from the FieldTrip-toolbox:
    % figure; hold on;
    % ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facecolor', 'skin'); alpha 0.5 
    % ft_plot_sens(sensors, 'elecshape', 'disc', 'elecsize', 1, 'label', 'off');
    % view([1 1 1])
    % axis on; grid on
    % xlabel('x'); ylabel('y'); zlabel('z'); 
    %  
    % ni2_topoplot(sensors, leadfield1); colorbar
    
    % % Just a check:
    % leadfield_OPM = mean(leadfield1(1:301));
    % leadfield_SQUID = mean(leadfield); 
    % 
    % % Conclusion: leadfield_SQUID > leadfield_OPM. I would expect the opposite
    % % since OPM are closer to the scalp ?? Probably is because of the
    % % different sensor noises
    
    leadfield_pinv_OPM_default=pinv(leadfield1(1:301)); % I take only the radial channels
    leadfield_pinv_OPM_default_mean = mean(abs(leadfield_pinv_OPM_default));
    
    sensor_number_OPM(i) = leadfield_pinv_OPM(i) / leadfield_pinv_OPM_default_mean;
end

% Plot the results
figure;
scatter(sensor_number_OPM, sensor_noise_OPM) % I expect sensor_noise_OPM ~ sqrt(sensor_number_OPM) but this does not happen. Why?
xlabel('OPM sensor number')
ylabel('OPM sensor noise')

figure;
scatter(sensor_number_OPM, leadfield_pinv_OPM) % This should be the same as the fig. 2E from the paper https://doi.org/10.1016/j.neuroimage.2021.118025. But it is not the same. Why?
xlabel('OPM sensor number')
ylabel('||pinv(L)||')

% Conclusion: I do NOT get the plots that I am expecting!




%% Way2: what if I use the leadfields and not the pinv

sensor_noise_OPM = 10 * 10^-10; % Fieldline paper: the remaining average noise floor (after SSP, look figure 5B) is 20 ± 5 fT/Hz^1/2.
sensor_noise_SQUID =  1 * 10^-10; % In intranet it says: The noise levels of MEG sensors should be below 10 fT, usually around 4-7fT. (https://intranet.donders.ru.nl/index.php?id=meg-lab-procedure-noiserecording&no_cache=1&sword_list%5B%5D=lab)

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors = ni2_sensors('type', 'ctf275');  
dippar = [0 0 6 1 0 0];
leadfield = ni2_leadfield(sensors, headmodel, dippar, sensor_noise_SQUID);

% Idea from Brookes et al. (2021) which is the same as eq. (4) (see photo in https://github.com/contsili/master-thesis/issues/2#issuecomment-1799608261)
leadfield_OPM = norm(leadfield) * sensor_noise_OPM / sensor_noise_SQUID;

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); % headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3) doesn't work since 'cocentric spheres' is not supported as a head model.
sensors = ni2_sensors('type', 'opm_tangential_radial');
dippar1 = [0 0 6 1 0 0]; 
leadfield1 = ni2_leadfield(sensors, headmodel, dippar1, sensor_noise_OPM);
leadfield_SQUID_default_mean = mean(abs(leadfield1(1:301)));

sensor_number_OPM = leadfield_OPM / leadfield_SQUID_default_mean







%% 14/11/2023
%% Sanity checks 1.,2.,3.

% sensor_noise_OPM = [1e-10, 5e-10, 7e-10, 9e-10, 1e-9, 1.1e-9, 1.3e-9, 1.6e-9, 2e-9, 2.5e-9]%, 4e-9];
% sensor_noise_OPM = linspace(0.001, 101, 20)*10^-10;
sensor_noise_OPM = 3*10^-10;

leadfield_OPM_eq4 = zeros(1,length(sensor_noise_OPM));

% Loop over different OPM sensor noise levels
for i = 1:length(sensor_noise_OPM)

    % Set a fixed SQUID sensor noise level based on lab information
    sensor_noise_SQUID =  1e-10;
 
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'ctf275');      
    dippar1 = [0 0 6 1 0 0];  
    leadfield = ni2_leadfield(sensors, headmodel, dippar1)
%     leadfield = leadfield + sensor_noise_SQUID; % Robert's suggestion: the noise should be stable and not change based on randn().


    % I want the dipole moment uncertainty of the OPM to be equal to the one from the SQUID. I compute that analytically using eq. (4) from https://github.com/contsili/master-thesis/issues/2#issuecomment-1799608261
    leadfield_OPM_eq4(i) = norm(leadfield) * sensor_noise_OPM(i) / sensor_noise_SQUID
    
    % Solve the forward model for a radial OPM system
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'opm_radial');    
    dippar2 = dippar1;     
    leadfield1 = ni2_leadfield(sensors, headmodel, dippar2); % Robert's suggestion: there should be no noise that changes based on randn(). So do not use: leadfield1 = ni2_leadfield(sensors, headmodel, dippar2, sensor_noise_OPM(i));
    leadfield_OPM_forward = norm(leadfield1)

    % Find the number of OPM sensors by weighting the norm of
    % leadfield from eq. 4 with the norm of the leadfield by forward model for a radial OPM system    
    sensor_number_OPM(i) = (sqrt(length(leadfield1)) * leadfield_OPM_eq4(i) / leadfield_OPM_forward)^2

    %% (Optional) The following code snippet was an attempt to calculate sensor number differently, but it does not give the expected results:
    % sensitivity = mean(leadfield1(1:301));
    % sensor_number_OPM(i) = (leadfield_OPM(i) / sensitivity)^2
    
end

% Plot the results
figure;
scatter(sensor_number_OPM, sensor_noise_OPM, 'r') % I expect sensor_noise_OPM ~ sqrt(sensor_number_OPM) but this does not happen. Why?
xlabel('OPM sensor number')
ylabel('OPM sensor noise')

% figure;
% scatter(sensor_number_OPM, leadfield_OPM_eq4) % This should be the same as the fig. 2E from the paper https://doi.org/10.1016/j.neuroimage.2021.118025. But it is not the same. Why?
% xlabel('OPM sensor number')
% ylabel('||L||')

% Some of my results are saved in
% C:\Users\user\Documents\Courses\Internship\master-thesis\simulations\results\aim2\analytical.
% These results are also in https://github.com/contsili/master-thesis/issues/2#issuecomment-1811313055





%% Sanity check 4.
sensor_noise_SQUID_unknown = 1 * 10^-15;
sensor_noise_SQUID =  1 * 10^-15; 

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors = ni2_sensors('type', 'ctf275');  
dippar = [0 0 6 1 0 0];
leadfield = ni2_leadfield(sensors, headmodel, dippar);

% I want the dipole moment uncertainty of the OPM to be equal to the one from the SQUID. I compute that analytically using eq. (4) from https://github.com/contsili/master-thesis/issues/2#issuecomment-1799608261    
leadfield_SQUID = norm(leadfield) * sensor_noise_SQUID_unknown / sensor_noise_SQUID;

% Compute 50 instances of the number of sensors for the SQUID system with
% an unknown number of sensors and in the average:
for i=1:50 

    % Solve the forward model for a SQUID system with an unknown number of
    % sensors
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'ctf275');      
    dippar = [0 0 6 1 0 0];   
    leadfield1 = ni2_leadfield(sensors, headmodel, dippar);
    leadfield_SQUID_default_mean(i) = norm(leadfield1);
    
    % Find the number of OPM sensors by weighting the norm of
    % leadfield from eq. 4 with the norm of the leadfield by forward model for a radial OPM system    
    sensor_number_SQUID(i) = (length(leadfield1)) * leadfield_SQUID / leadfield_SQUID_default_mean(i)
end

sensor_number_SQUID_average = mean(sensor_number_SQUID)

% 20 repetitions
% 1st run: sensor_number_OPM_average = 275.6220
% 2nd: sensor_number_OPM_average = 281.4209
% 3rd: sensor_number_OPM_average = 262.2722

% 100 repetitions
% 4th: sensor_number_OPM_average = 269.5290
% 5th: sensor_number_OPM_average = 282.0048
% 6th: sensor_number_OPM_average = 273.3936



%% An alternative to "sensor_number_OPM(i) = length(leadfield1) * leadfield_OPM_eq4(i) / leadfield_OPM_forward": explicitly check norm(leadfield_OPM_eq4) == norm(leadield_OPM_forward)

sensor_noise_OPM = [5 * 10^-10, 7 * 10^-10, 9 * 10^-10, 10 * 10^-10, 11 * 10^-10, 13 * 10^-10, 16 * 10^-10,  20 * 10^-10, 25 * 10^-10];

% SQUID
for i = 1:length(sensor_noise_OPM)
    sensor_noise_SQUID =  1 * 10^-10; % In intranet it says: The noise levels of MEG sensors should be below 10 fT, usually around 4-7fT. (https://intranet.donders.ru.nl/index.php?id=meg-lab-procedure-noiserecording&no_cache=1&sword_list%5B%5D=lab)
      
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'ctf275');   
    dippar = [0 0 6 1 0 0];    
    leadfield = ni2_leadfield(sensors, headmodel, dippar, sensor_noise_SQUID);
        
    % Idea from Brookes et al. (2021) which is the same as eq. (4) (see photo in https://github.com/contsili/master-thesis/issues/2#issuecomment-1799608261)    
    leadfield_OPM(i) = norm(leadfield) * sensor_noise_OPM(i) / sensor_noise_SQUID;
    
    %% 1. opm_radial
    
    sensor_number_OPM = 21:10:350;

    for j=1:length(sensor_number_OPM)
        headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); % headmodel = ni2_headmodel('type', 'spherical', 'nshell', 3) doesn't work since 'cocentric spheres' is not supported as a head model.
        sensors = ni2_sensors('type', 'opm_tangential_radial','sensor_number', sensor_number_OPM(j));
        
        dippar1 = [0 0 6 1 0 0]; 
        
        leadfield1 = ni2_leadfield(sensors, headmodel, dippar1);

        if round(norm(leadfield1(1:length(leadfield1)/2)),8)==round(leadfield_OPM(i),8)

            final_sensor_number=sensor_number_OPM(j)
            break;
        end
    end
end

% IT DOES NOT WORK





%% 17/11/2023 (just before the OPM installation) 

% 1. Here I use Robert's suggestion: ||pinv(L)|| ~= 1/||L||. So, in eq. (4) use ||pinv(L)|| (not ||L||) - In that case do not trust the Brooke's paper
% Conclusion to 1.: the results with ||pinv(L)|| (17/11/2023) were as bad as with ||L|| (14/11/2023).

% 2. Since 1. does not work I try a different trick which is "methodos twn triwn" with norm(pinv(leadfield1) ~ 1/ sqrt(number of sensors)



%% Sanity checks 1.,2.,3. 
% sensor_noise_OPM = [1e-10, 5e-10, 7e-10, 9e-10, 1e-9, 1.1e-9, 1.3e-9, 1.6e-9, 2e-9, 2.5e-9]%, 4e-9];
% sensor_noise_OPM = linspace(0.001, 100, 20)*10^-10;
sensor_noise_OPM = 3e-10;
leadfield_OPM_eq4 = zeros(1,length(sensor_noise_OPM));

% Set a fixed SQUID sensor noise level based on lab information
sensor_noise_SQUID =  1e-10;


headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors = ni2_sensors('type', 'ctf275');      
dippar1 = [0 0 6 1 0 0];  
leadfield = ni2_leadfield(sensors, headmodel, dippar1);
pinv_leadfield = pinv(leadfield);

% Loop over different OPM sensor noise levels
for i = 1:length(sensor_noise_OPM)

    
    % I want the dipole moment uncertainty of the OPM to be equal to the one from the SQUID. I compute that analytically using eq. (4) from https://github.com/contsili/master-thesis/issues/2#issuecomment-1799608261
    norm_pinv_leadfield_OPM_eq4(i) = norm(pinv_leadfield) * sensor_noise_SQUID / sensor_noise_OPM(i);
    
    % Both ||pinv(L)|| and ||L|| should depend on sqrt(number of OPM
    % sensors). Note: ||pinv(L)|| ~ sqrt(number of OPM sensors) can not
    % hold because then eq. 4 gives that more noise leads to less sensors
    % (which cannot hold)

    % Solve the forward model for a radial OPM system
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'opm_radial');    
    dippar2 = dippar1;     
    leadfield1 = ni2_leadfield(sensors, headmodel, dippar2);
    norm_pinv_leadfield_OPM_forward = norm(pinv(leadfield1));

    % Find the number of OPM sensors by weighting the norm of
    % leadfield from eq. 4 with the norm of the leadfield by forward model for a radial OPM system    
    sensor_number_OPM(i) = (sqrt(length(leadfield1)) * norm_pinv_leadfield_OPM_forward / norm_pinv_leadfield_OPM_eq4(i) )^2 

    %% Run some checks:

    %% Check1
    % The equation in the previous line holds since: ||pinv(L)|| ~ 1/sqrt(number of OPM sensors). Lets do a check:
    %    >> norm_pinv_leadfield_OPM_forward = norm(pinv(leadfield1(1:200)))
    % 
    % norm_pinv_leadfield_OPM_forward =
    % 
    %    8.2965e+07
    % 
    %    >> norm_pinv_leadfield_OPM_forward = norm(pinv(leadfield1(1:100)))
    % 
    % norm_pinv_leadfield_OPM_forward =
    % 
    %    1.3072e+08
    %  
    %   Let's say that we do not know  the number of sensors that the
    %   norm_pinv_leadfield_OPM_forward=1.3072e+08 has:
    %   8.2965e+07 -> 1/sqrt(200)
    %   1.3072e+08 -> ? = 0.1114 (note 1/sqrt(100) = 0.1000. It is not totally the same since leadfield1 does not have the same value for every element/sensor)

    %% Check2
    %  Also the plot I create below looks like it changes 1/sqrt(number of sensors):
    % for i=10:10:200
    % scatter(i,norm(pinv(leadfield1(1:i))))
    % hold on 
    % end
    
    %% Check3
%     % Question: Can you check if the relation between i and norm(pinv(leadfield1(1:i)) is 1/x or 1/sqrt(x)
%     % Answer: let matlab search for the best solution for n when y=1/x^n
%     x = [10:10:200]';
%     
%     for i = 1:length(x)
%         y(i) =  norm(pinv(leadfield1(1:x(i))));
%     end
%     y=y';
%     
%     b=polyfit(log10(x),log10(y),1);
%     b(1)=-b(1); b(2)=10^b(2);   % for initial start point in nonlinear model
%     fo=fitoptions('Method','NonlinearLeastSquares','StartPoint',flip(b));
%     ft=fittype(@(a,n,x) a./x.^n,'options',fo);
%     
%     f1=fit(x,y,ft);
%     
%     %  As a result I get that y=a/x^n, where n = 0.5416  (CI is: [0.5029,
%     %  0.5803]) which is close to 1/2. So  norm(pinv(leadfield1)) =
%     %  a/x^(1/2)=a/sqrt(x)
    
end

% Plot the results
figure;
scatter(sensor_number_OPM, sensor_noise_OPM) % I expect sensor_noise_OPM ~ sqrt(sensor_number_OPM) but this does not happen. Why?
xlabel('OPM sensor number')
ylabel('OPM sensor noise')

figure;
scatter(sensor_number_OPM, norm_pinv_leadfield_OPM_eq4) % This should be the same as the fig. 2E from the paper https://doi.org/10.1016/j.neuroimage.2021.118025. But it is not the same. Why?
xlabel('OPM sensor number')
ylabel('||pinv(L)||')


% Conclusion: IT WORKS! The grahp is ||L|| = sqrt(sensor_number_OPM)

%% 08/01/2023

% 1. The relationship between ||L|| and sqrt(sensor_number_OPM) is not necessarily an sqrt(). It is quite random since it depends on what sensor I add each time. ||L|| is high when the sensor close to the source and low when the sensor is far from the source.
% Solution to 1.: I use non-linear least square fitting (17/11/2023) to find the ideal n and I do not force ||pinv(L)|| = 1/sqrt(sensor_number_OPM)

%% Sanity checks 1.,2.,3. 
sensor_noise_OPM = [1e-10, 5e-10, 7e-10, 9e-10, 1e-9, 1.1e-9, 1.3e-9, 1.6e-9, 2e-9, 2.5e-9];%, 4e-9];
% sensor_noise_OPM = linspace(0.001, 100, 20)*10^-10;
% sensor_noise_OPM = 3e-10;
leadfield_OPM_eq4 = zeros(1,length(sensor_noise_OPM));

% Set a fixed SQUID sensor noise level based on lab information
sensor_noise_SQUID =  1e-10;


headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors = ni2_sensors('type', 'ctf275');      
dippar1 = [0 0 6 1 0 0];  
leadfield = ni2_leadfield(sensors.ctf275, headmodel, dippar1);
pinv_leadfield = pinv(leadfield);

% Loop over different OPM sensor noise levels
for i = 1:length(sensor_noise_OPM)

    norm_pinv_leadfield_OPM_eq4(i) = norm(pinv_leadfield) * sensor_noise_SQUID / sensor_noise_OPM(i);

    % Solve the forward model for a radial OPM system
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'opm_radial');    
    dippar2 = dippar1;     
    leadfield1 = ni2_leadfield(sensors, headmodel, dippar2);
    norm_pinv_leadfield_OPM_forward = norm(pinv(leadfield1));
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Question: Can you check if the relation between i and norm(pinv(leadfield1(1:i)) is 1/x or 1/sqrt(x)
    % Answer: let matlab search for the best solution for n when
    % norm(pinv(leadfield1))=1/sensor_number_OPM^n. In other words, do
    % nonlinear least squares fitting. The goal is to fit a model function to the data in such a way that it minimizes the sum of the squares of the differences between the observed (number of sensors) and predicted values (norm(pinv(leadfield))).
     
    % Find the index of the first non-zero element in leadfield1 
    k = find(leadfield1, 1, 'first');
    
    % Create the vector x starting from the index k
    x = (k:length(leadfield1))'; % if I put x = (1:length(leadfield1))'; then I get wrong calculation of n since leadfield1(1) might be a equal to 0
    
    y = zeros(length(x),1);
    for l = 1:length(x)
        y(l) =  norm(pinv(leadfield1(1:x(l))));
    end
    y=y(:);
    
    b=polyfit(log10(x),log10(y),1);
    b(1)=-b(1); b(2)=10^b(2);   % for initial start point in nonlinear model
    fo=fitoptions('Method','NonlinearLeastSquares','StartPoint',flip(b));
    ft=fittype(@(a,n,x) a./x.^n,'options',fo);
    
    f1=fit(x,y,ft);
    
    % Access the coefficients
    coefficients = coeffvalues(f1);
    
    % Extract the value of n
    n_value(i) = coefficients(2);
    %%%%%%%%%%%%%%%%%%%%%

    % Find the number of OPM sensors by weighting the norm of
    % leadfield from eq. 4 with the norm of the leadfield by forward model for a radial OPM system    
    sensor_number_OPM(i) = ( (length(leadfield1)^n_value(i)) * norm_pinv_leadfield_OPM_forward / norm_pinv_leadfield_OPM_eq4(i) )^(1/n_value(i)); 
end

% Plot the results
figure;
scatter(sensor_number_OPM, sensor_noise_OPM) % I expect sensor_noise_OPM ~ sqrt(sensor_number_OPM) but this does not happen. Why?
xlabel('OPM sensor number')
ylabel('OPM sensor noise')

figure;
scatter(sensor_number_OPM, norm_pinv_leadfield_OPM_eq4) % This should be the same as the fig. 2E from the paper https://doi.org/10.1016/j.neuroimage.2021.118025. But it is not the same. Why?
xlabel('OPM sensor number')
ylabel('||pinv(L)||')

% Conclusion: IT WORKS PERFECTLY!!!

% Note: I checked and n_value is the same in every loop!!


%% Do what I did before but use only the ||L|| not ||pinv(L)||

sensor_noise_OPM = [1e-10, 5e-10, 7e-10, 9e-10, 1e-9, 1.1e-9, 1.3e-9, 1.6e-9, 2e-9, 2.5e-9];%, 4e-9];
% sensor_noise_OPM = linspace(0.001, 100, 20)*10^-10;
% sensor_noise_OPM = 3e-10;
leadfield_OPM_eq4 = zeros(1,length(sensor_noise_OPM));

% Set a fixed SQUID sensor noise level based on lab information
sensor_noise_SQUID =  1e-10;

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors = ni2_sensors('type', 'ctf275');      
dippar1 = [0 0 6 1 0 0];  
leadfield = ni2_leadfield(sensors.ctf275, headmodel, dippar1);
pinv_leadfield = pinv(leadfield);

for i = 1:length(sensor_noise_OPM)
    
    norm_leadfield_OPM_eq4(i) =  sensor_noise_OPM(i) *  norm(leadfield)/  sensor_noise_SQUID; 
    
    % Solve the forward model for a radial OPM system
    headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
    sensors = ni2_sensors('type', 'opm_radial');    
    dippar2 = dippar1;     
    leadfield1 = ni2_leadfield(sensors, headmodel, dippar2);
    norm_leadfield_OPM_forward = norm(leadfield1);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Question: Can you check if the relation between i and norm(pinv(leadfield1(1:i)) is 1/x or 1/sqrt(x)
    % Answer: let matlab search for the best solution for n when
    % norm(pinv(leadfield1))=1/sensor_number_OPM^n. In other words, do
    % nonlinear least squares fitting. The goal is to fit a model function to the data in such a way that it minimizes the sum of the squares of the differences between the observed (number of sensors) and predicted values (norm(pinv(leadfield))).
     
    % Find the index of the first non-zero element in leadfield1 
    k = find(leadfield1, 1, 'first');
    
    % Create the vector x starting from the index k
    x = (k:length(leadfield1))'; % if I put x = (1:length(leadfield1))'; then I get wrong calculation of n since leadfield1(1) might be a equal to 0
    x=x(:);

    y = zeros(length(x),1);
    for l = 1:length(x)
        y(l) =  norm(leadfield1(1:x(l)));
    end
    y=y(:);
       
    power_law_model = fittype('a * x^n', 'coefficients', {'a', 'n'});
    fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [mean(y), 0.5]);
    fitted_model = fit(x, y, power_law_model, fit_options);
    
    % Access the fitted coefficients
    a = fitted_model.a;
    n = fitted_model.n;
    %%%%%%%%%%%%%%%%%%%%%

    % way1
    sensor_number_OPM(i) = ( (length(leadfield1)^n) * norm_leadfield_OPM_eq4(i) / norm_leadfield_OPM_forward )^(1/n); 
    
    % way2   
%     % Find the value in y that is closest to sensor_number_OPM. The sensor_number_OPM cannot
%     % go more than 299
%     [~, index] = min(abs(y - norm_leadfield_OPM_eq4(i)));
%     sensor_number_OPM(i) = x(index);
%     
%     fprintf('Closest value: %f\n', sensor_number_OPM);
%     fprintf('Index of closest value: %d\n', index);
end

figure;
scatter(sensor_number_OPM, sensor_noise_OPM) % I expect sensor_noise_OPM ~ sqrt(sensor_number_OPM) but this does not happen. Why?
xlabel('OPM sensor number')
ylabel('OPM sensor noise')
hold off

figure;
y1=a*x.^(n);
plot(x, y, 'o', 'DisplayName', 'Data');
hold on;
plot(x, y1, 'r-', 'DisplayName', 'Fit');

