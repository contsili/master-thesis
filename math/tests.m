
%% TEST 1
% What if I take less sensors for the forward model? 
% Expected result: I should get the same sensor_number_OPM as when I have the full 301 sensors forward model. In other words the forward model should not affect the relationship between ||pinv(L)|| and 1/sqrt(sensor_number_OPM) 

% 1. I should use mesh_sphere() (inside the ni2_sensors()) for that,
% otherwise in each loop I will get a different coverage of the sensors. So
% I use  sensors = ni2_sensors('type', 'eeg', 'n', 100);

% 09/01/2024: THIS TEST IS NOT NEEDED! I am just fitting a line to find the
% best approximation for
% n_value, so the more points the better (so keep the max points which are
% 301 sensors)

% Meeting with Robert 09/01/2024: To have ||pinv(L)|| depend only on number of sensors: 
% 1. for the forward opm model I need to keep the sensors in a stable helmet (I do that)
% 2. I need to keep the dipole in a certain position (I do that)
% 3. I need to not use noise in the ni2_leadfield. This randn component will make my model unnecessarily complex (I do that)


%% TEST 2
% OPM: Put the OPMs on a SQUID helmet. Give the OPMs sensor_noise_OPM = sensor_noise_SQUID
% SQUID: For the SQUID use only the magnetometers 
% Expected result: sensor_number_OPM = 273 in order for OPM to have the
% same performance with SQUID


%% TEST 3 
% Do not force eq.4 to be equal to 1 AND keep the same noise for SQUID and
% OPMs (N_opm=N_squid). Similar to method 2a from analytical_RO.m 

% ctf275 vs opm301

% SQUID (has 275 gradiometers channels)
headmodel   = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors     = ni2_sensors('type', 'ctf275'); % Alternative: load C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\ctf275.mat and then only keep ctf275.chantype = meggrad.     
dippar      = [0 0 6 1 0 0];  
leadfield   = ni2_leadfield(sensors.ctf275, headmodel, dippar);

L_squid     = norm(leadfield);
pinvL_squid = norm(pinv(leadfield));

% OPM (has 301 radial magnetometer channels)
headmodel  = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors    = ni2_sensors('type', 'opm_radial');    
dippar     = [0 0 6 1 0 0];  
leadfield  = ni2_leadfield(sensors, headmodel, dippar);

L_opm      = norm(leadfield); 
pinvL_opm  = norm(pinv(leadfield));

% eq. 4 (N_opm=N_squid, so they cancel each other in the equation)

L_opm / L_squid 

L_opm / (3*L_squid) 

% Conclusion: opm better than squid as expected
% 
% L_opm / L_squid = 5.0454. If we add sensor noise (N_opm=3*N_squid) then 
% (L_opm / N_opm) / (L_squid / N_squid) =  L_opm / 3 * L_squid. 
% 
% Also if we consider the 144 channels from the fiedlinebeta2 we have
% approximately half channels than the 301 of the current helmet, so: 
% ( L_opm) / (3* sqrt(2) * L_squid) ~ 1.18. 


%% TEST 4a: ctf275 vs fieldlinebeta2 for DEEP dipole -> OPM loses

sensor_noise_SQUID =  1e-10;
sensor_noise_OPM   = 3 * sensor_noise_SQUID;
ntime              = 1000;
s                  = ones(1,ntime); % random choise of the time course of source
dippar             = [0 0 6 1 0 0];  

% SQUID (has 275 gradiometers channels)
headmodel   = ni2_headmodel('type', 'spherical', 'nshell', 1); 
load C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\ctf275.mat 
sel275      = startsWith(ctf275.label, 'M');   
lf275       = ni2_leadfield(ctf275, headmodel, dippar);
lf275       = lf275(sel275, :);
norm_lf275  = norm(lf275);

dat275 = lf275 * s + sensor_noise_SQUID * randn(275,ntime);
est275 = pinv(lf275) * dat275;

% fieldlinebeta2 with 432 channels
headmodel  = ni2_headmodel('type', 'spherical', 'nshell', 1); 
load C:\Users\user\Documents\Courses\Internship\master-thesis\simulations\ni2-electrophys-master\fieldlinebeta2\grad.mat;  
lf432      = ni2_leadfield(grad, headmodel, dippar);
norm_lf432 = norm(lf432); 

dat432 = lf432 * s + sensor_noise_OPM * randn(432,ntime);
est432 = pinv(lf432) * dat432;

% fieldlinebeta2 with 144 channels
headmodel  = ni2_headmodel('type', 'spherical', 'nshell', 1); 
load C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\fieldlinebeta2.mat
lf144      = ni2_leadfield(fieldlinebeta2, headmodel, dippar);
norm_lf144 = norm(lf144); 

dat144 = lf144 * s + sensor_noise_OPM * randn(144,ntime);
est144 = pinv(lf144) * dat144;

% ctf275 vs fl432
disp('fl432:')
norm_lf432 / norm_lf275 % Wihtout sensor noise: 2.3091
std(est275) / std(est432) % With sensor noise: 0.7801

% ctf275 vs fl144
disp('fl144:')
norm_lf144 / norm_lf275 % % Without sensor noise: 1.7941
std(est275) / std(est144) % % With sensor noise: 0.6176  
% L_fl144 / (3*L_squid) % this is equal to std(est275) / std(est144) because of eq. 4

% plot the sensors, the head model and the dipole
figure; 
ft_plot_sens(ctf275, 'label', 'no', 'axes', 0, 'orientation', 0, 'marker','.','unit', 'cm');
hold on
ft_plot_headmodel(headmodel,'facealpha', 0.7)
ft_plot_dipole(dippar(1:3), dippar(4:6))

figure;
ft_plot_sens(grad, 'label', 'no', 'axes', 0, 'orientation', 0, 'marker','o','unit', 'cm', 'edgecolor', 'k');
hold on 
ft_plot_headmodel(headmodel, 'facealpha', 0.7)
ft_plot_dipole(dippar(1:3), dippar(4:6))

figure
ft_plot_topo3d(fieldlinebeta2.chanpos, lf144)
figure
ft_plot_topo3d(grad.chanpos, lf432)
figure
ft_plot_topo3d(ctf275.chanpos(sel275,:), lf275)

%% TEST 4b: ctf275 vs fieldlinebeta2 for SUPERFICIAL dipole -> OPM wins 

sensor_noise_SQUID =  1e-10;
sensor_noise_OPM   = 3 * sensor_noise_SQUID;
ntime              = 1000;
s                  = ones(1,ntime); % random choise of the time course of source
dippar             = [0 0 10 1 0 0];  

% SQUID (has 275 gradiometers channels)
headmodel   = ni2_headmodel('type', 'spherical', 'nshell', 1); 
load C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\ctf275.mat 
sel275      = startsWith(ctf275.label, 'M');   
lf275       = ni2_leadfield(ctf275, headmodel, dippar);
lf275       = lf275(sel275, :);
norm_lf275  = norm(lf275);

dat275 = lf275 * s + sensor_noise_SQUID * randn(275,ntime);
est275 = pinv(lf275) * dat275;

% fieldlinebeta2 with 432 channels
headmodel  = ni2_headmodel('type', 'spherical', 'nshell', 1); 
load C:\Users\user\Documents\Courses\Internship\master-thesis\simulations\ni2-electrophys-master\fieldlinebeta2\grad.mat;  
lf432      = ni2_leadfield(grad, headmodel, dippar);
norm_lf432 = norm(lf432); 

dat432 = lf432 * s + sensor_noise_OPM * randn(432,ntime);
est432 = pinv(lf432) * dat432;

% fieldlinebeta2 with 144 channels
headmodel  = ni2_headmodel('type', 'spherical', 'nshell', 1); 
load C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\fieldlinebeta2.mat
lf144      = ni2_leadfield(fieldlinebeta2, headmodel, dippar);
norm_lf144 = norm(lf144); 

dat144 = lf144 * s + sensor_noise_OPM * randn(144,ntime);
est144 = pinv(lf144) * dat144;

% ctf275 vs fl432
disp('fl432:')
norm_lf432 / norm_lf275 % Wihtout sensor noise: 2.3091
std(est275) / std(est432) % With sensor noise: 0.7801

% ctf275 vs fl144
disp('fl144:')
norm_lf144 / norm_lf275 % % Without sensor noise: 1.7941
std(est275) / std(est144) % % With sensor noise: 0.6176  
% L_fl144 / (3*L_squid) % this is equal to std(est275) / std(est144) because of eq. 4

% plot the sensors, the head model and the dipole
figure; 
ft_plot_sens(ctf275, 'label', 'no', 'axes', 0, 'orientation', 0, 'marker','.','unit', 'cm');
hold on
ft_plot_headmodel(headmodel,'facealpha', 0.7)
ft_plot_dipole(dippar(1:3), dippar(4:6))

figure;
ft_plot_sens(grad, 'label', 'no', 'axes', 0, 'orientation', 0, 'marker','o','unit', 'cm', 'edgecolor', 'k');
hold on 
ft_plot_headmodel(headmodel, 'facealpha', 0.7)
ft_plot_dipole(dippar(1:3), dippar(4:6))

figure
ft_plot_topo3d(fieldlinebeta2.chanpos, lf144)
figure
ft_plot_topo3d(grad.chanpos, lf432)
figure
ft_plot_topo3d(ctf275.chanpos(sel275,:), lf275)

%% TEST 5
% plot ||L|| - number of sensors for the fieldlinebeta2. It looks like a "staircase"

headmodel  = ni2_headmodel('type', 'spherical', 'nshell', 1); 
load C:\Users\user\Documents\Courses\Internship\master-thesis\simulations\ni2-electrophys-master\fieldlinebeta2\grad.mat;  
dippar      = [0 0 6 0 1 0];  
leadfield   = ni2_leadfield(grad, headmodel, dippar);

% Non-linear least square fitting
k = find(leadfield, 1, 'first');

x = (k:length(leadfield))';
x=x(:);

y = zeros(length(x),1);
for l = 1:length(x)
    y(l) =  norm(leadfield(1:x(l)));
end
y=y(:);

power_law_model = fittype('a * x^n', 'coefficients', {'a', 'n'});
fit_options = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [y(1), 0.52]);
fitted_model = fit(x, y, power_law_model, fit_options);

a = fitted_model.a;
n = fitted_model.n;

% Plot
figure;
y1=a*x.^(n);
plot(x, y, 'o', 'DisplayName', 'Data');
hold on;
plot(x, y1, 'r-', 'DisplayName', 'Fit');
xlabel('OPM sensor number')
ylabel('||L||')