%% Simulated data

%% JK
% Description: I fit a dipole for 146 trials.

%% this is the code I run in the local PC

addpath('H:/common/matlab/fieldtrip')

ft_defaults

%% OPM

load M:\Documents\recordings\opm-MN\results\coregistration\hfc\'mat files'\headmodel_sphere.mat
load M:\Documents\recordings\opm-MN\results\dipmom_analytical\mom135.mat
load M:\Documents\recordings\opm-MN\results\dipmom_analytical\pos135.mat
load M:\Documents\recordings\opm-MN\results\jackknife\hfc\fieldlinebeta2_head_fixedsensors.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\dpu_squid.mat
load M:\Documents\recordings\opm-MN\results\jackknife\hfc\data_meg_hfc_rejectvisual_146trials.mat

load M:\Documents\recordings\opm-MN\results\dippos_jk\squid\squid_noise_floor.mat


%% Prepare the data

% % To be able to do jackknifing I need all the runs to have the same number of trials. I will keep 176 trials for all the runs
% 
% for i=1:6
%     cfg        = [];
%     cfg.method = 'summary';
%     data_meg_hfc_clean(i) = ft_rejectvisual(cfg, data_meg_hfc(i)); 
% end



%% Rename sensor that re-appear between runs

% I had 8 stable sensors between runs. I will rename them
for i = 2:6     
    index = find(ismember(data_meg_hfc_clean(i).label, {'R101_bz', 'L108_bz','R503_bz','L503_bz','L507_bz','R507_bz','R212_bz','L212_bz'}));
    for k = 1:length(index)
        data_meg_hfc_clean(i).label{index(k)} = strcat(data_meg_hfc_clean(i).label{index(k)}, '_', num2str(i));
    end    
end

% the 6th run has some channels that appeared in the 5th run too. I am going to
% rename it.
index = find(ismember( data_meg_hfc_clean(6).label,  data_meg_hfc_clean(5).label));

for k = 1:length(index)   
        data_meg_hfc_clean(6).label{index(k)} = strcat(data_meg_hfc_clean(6).label{index(k)}, '-', num2str(6));
end

% Also R114 appears in the 6th and 1st run so I rename it.
index_R114 = find(strcmp(data_meg_hfc_clean(6).label, 'R114_bz'));

if ~isempty(index_R114)
    data_meg_hfc_clean(6).label{index_R114} = 'R114_bz-6';
end


% After I checked all the channels that I know that they are duplicate, now I will see which channels still remain duplicate        
% Now I rename these duplicate channels with the logic: L214_bz is
% duplicate at Run 2 and Run 3. So, I rename the run 3 as L214_bz__3:


data_meg_hfc_clean(5).label{16} = strcat(data_meg_hfc_clean(5).label{16}, '__5');
data_meg_hfc_clean(5).label{31} = strcat(data_meg_hfc_clean(5).label{31}, '__5');
data_meg_hfc_clean(3).label{27} = strcat(data_meg_hfc_clean(3).label{27}, '__3');


%% After I checked all the channels that I know that they are duplicate, now I will see which channels still remain duplicate
duplicate_channels = {};

for i = 1:6
    % Loop through each channel in the current run
    for j = 1:length(data_meg_hfc_clean(i).label)
        current_channel = data_meg_hfc_clean(i).label{j};
        % Compare the current channel with channels in other runs
        for k = (i + 1):6
            % Check if the channel name exists in the other runs
            if any(strcmp(current_channel, data_meg_hfc_clean(k).label))
                % If found duplicate, add it to the duplicate_channels cell array
                duplicate_channels{end + 1} = {current_channel, ['Run ', num2str(i), ' and Run ', num2str(k)]};
            end
        end
    end
end

% Display the duplicate channels
if ~isempty(duplicate_channels)
    disp('Duplicate Channels:');
    disp('Channel Name       :   Runs');
    for m = 1:length(duplicate_channels)
        disp([duplicate_channels{m}{1}, '   :   ', duplicate_channels{m}{2}]);
    end
else
    disp('No duplicate channels found.');
end


%% append data

cfg         = [];
cfg.keepsampleinfo='no';
append_data = ft_appenddata(cfg, data_meg_hfc_clean(1), data_meg_hfc_clean(2), data_meg_hfc_clean(3), data_meg_hfc_clean(4), data_meg_hfc_clean(5), data_meg_hfc_clean(6)); 

append_data.grad = fieldlinebeta2_head_fixedsensors; % QUESTION: Is is correct? Do I need to create a grad structure with 182 sensors as my sensors in append_data?

% NOTE: I have not deleted  cfg.channel     = {'all', '-L101_bz',
% '-L214_bz'}; that were not coregistered well. I think that
% fieldlinebeta2_head_fixedsensors will take care of that

%% I need to make 176 append_tl_jk.trial = L*q + add stochastic noise on my own. So basically get a 1 q and then ft_prepare_leadfield 176 times
% If I find that the relation of dpu and nr_sens is 1/sqrt() and see how it
% depends on the noise then we could use an analytical solution for our
% whole brain simulations

n = length(data_meg_hfc_clean(i).trial);


% sensor noise
% sensor_noise_squid = 5 * 10^-15;

sensor_noise_opm =  [1e-11,1e-12,1e-13,1e-14,1e-15,1e-16];


% cfg            = [];
% cfg.channel    = {'all', '-L101_bz', '-L214_bz'}; % delete the channels that appeared closer to the brain during the co-registration
% append_tl  = ft_timelockanalysis(cfg, append_data);
% 
% cfg            = [];
% cfg.latency    = [0.035 0.05];
% cfg.unit       = 'cm';
% cfg.gridsearch = 'yes';
% cfg.nonlinear  = 'yes'; 
% cfg.grad       = fieldlinebeta2_head_fixedsensors;                      
% cfg.headmodel  = headmodel_sphere;
% stim_dpN20_opm = ft_dipolefitting(cfg, append_tl);

% mom135_146trials = mean(stim_dpN20_opm.dip.mom, 2);

%% compute leadfield
cfg                 = [];      
cfg.sourcemodel.pos = pos135;
cfg.grad            = fieldlinebeta2_head_fixedsensors;
cfg.channel         = fieldlinebeta2_head_fixedsensors.label(1:128);
cfg.headmodel       = headmodel_sphere;
leadfield           = ft_prepare_leadfield(cfg);


mom = mom135; % data.avg (simulated signal) needs lf. lf needs orient (dipole orientation). orient needs mom (dipole moment). In an ideal simulation I generate a dipole and I know its mom. Here however, the best estimation I have is mom135 (estimation of moment from 135 OPM sensors)
norm_mom = norm(mom);
orient = mom ./ norm_mom;


%%
for l=1:length(sensor_noise_opm)
    % generate simulated data with sensor_noise_opm = 3*sensor_noise_squid
    new_data       = append_data;
    new_data.label = fieldlinebeta2_head_fixedsensors.label(1:128);
    for k = 1:146    
        new_data.trial{1,k} = leadfield.leadfield{1,1} * mom135 + sensor_noise_opm(l) * randn(length(fieldlinebeta2_head_fixedsensors.label(1:128)), 0.6*5000); % note: here I add the leadfield that is defined for 35-50 ms (when the dipolar pattern is there) to the whole trial length ((-200 ms) - 400 ms). This is not correct since the dipolar pattern appears only at 35-50 ms. However I correct that later when I use cfg.latency = [0.035 0.05];
    end
    for j = 1:n
         
        cfg            = [];
        cfg.trials     = setdiff(1:n, j); % leave one trial out
        % cfg.channel    = {'all', '-L101_bz', '-L214_bz'}; % delete the channels that appeared closer to the brain during the co-registration
        stim_tl_jk(j)  = ft_timelockanalysis(cfg, new_data);
        
        cfg            = [];
%         cfg.latency    = [0.035 0.05];
        cfg.latency    = 0.035;
        cfg.unit       = 'cm';
        cfg.gridsearch = 'no';
        cfg.dip.pos    = pos135; % location where dipole was localised with 135 OPM sensors. This is the starting point for the non-linear search
        cfg.grad       = fieldlinebeta2_head_fixedsensors;
%         cfg.channel    = fieldlinebeta2_head_fixedsensors.label(1:128);
        cfg.headmodel  = headmodel_sphere;
        stim_dpN20_opm_jk(:,j) = ft_dipolefitting(cfg, stim_tl_jk(j));
        
        dippos(j,:) = stim_dpN20_opm_jk(:,j).dip.pos;  

        % dmu
        dipmom(l,j) = pinv(leadfield.leadfield{1,1}*orient)*stim_dpN20_opm_jk(:,j).Vdata;
    
        
    end
    
    bias = (n-1)^2;
        
    jackknife_std_x = sqrt(bias) * std(dippos(:,1));
    jackknife_std_y = sqrt(bias) * std(dippos(:,2));
    jackknife_std_z = sqrt(bias) * std(dippos(:,3));
            
    dpu1(l) = sqrt((jackknife_std_x^2 + jackknife_std_y^2 + jackknife_std_z^2)/3)


    % dmu
    dmu1(l) = sqrt(bias)*sqrt(sum((dipmom(l,:)-norm_mom).^2)/146);
    
    dmu1_no_inflation(l) = sqrt(bias)*sqrt(sum((dipmom(l,:)-mean(dipmom(l,:), 'all')).^2)/146);

     
end

%%
figure
plot(log10(sensor_noise_opm), log10(dpu_var),'o')
xlabel('log10(noise)')
ylabel('log10(dpu)')




%% Monte-Carlo
% Description: I fit a dipole for 1 trial. 

% I know that noise lowers by 1/sqrt(trials). So to make the two methods
% the same:
sensor_noise_opm2 = sensor_noise_opm; 
  
%% compute leadfield



for l=1:length(sensor_noise_opm)

%     %%%%% Each sensor has each own sensor noise that does not change when
%     %%%%% we add another sensor on top of it %%%%%%%
%     rng(1, 'twister') % produce the same set of random numbers in each iteration
    
    
    %%%%%%%%%%% Dipole fitting %%%%%%%%%%%%%
    data        = [];
    data.time   = linspace(0.01,1,100); 
    data.avg    = leadfield.leadfield{1,1} * mom135 + sensor_noise_opm2(l) * randn(length(1:128), length(data.time));   %%%%%%%% Add Gaussian noise to leadfield %%%%%%%
    data.label  = fieldlinebeta2_head_fixedsensors.label(1:128);
    data.grad   = fieldlinebeta2_head_fixedsensors; 
    data.dimord = 'chan_time';
    
    cfg = [];
    cfg.gridsearch = 'no';
    cfg.dip.pos    = pos135;
    cfg.model      = 'moving';
    cfg.latency    = 'all';
    cfg.headmodel  = headmodel_sphere;
    cfg.nonlinear  = 'yes';
    cfg.numdipoles = 1;
    cfg.showcallinfo = 'no';
    
    dip = ft_dipolefitting(cfg, data); 
    
   
    from_structure_to_vector = extractfield(dip.dip, 'pos');
    dip_location_opm = reshape(from_structure_to_vector, [3,length(data.time)])'; 
    
    dpu2(l) = sqrt((std(dip_location_opm(:,1))^2 + std(dip_location_opm(:,2))^2 + std(dip_location_opm(:,3))^2)  / 3 ); % this formula is taken from Vrba (2000). "Multichannel SQUID biomagnetic systems"   

    
    % dmu2
    dipmom2(l,:) = pinv(leadfield.leadfield{1,1}*orient)*dip.Vdata;

    dmu2(l) = sqrt(sum((dipmom2(l,:)-norm_mom).^2)/100);

    dmu2_no_inflation(l) = sqrt(sum((dipmom2(l,:)-mean(dipmom2(l,:), 'all')).^2)/100);


end 

%% math

sensor_noise_opm3 = sensor_noise_opm; 


for l = 1:length(sensor_noise_opm)    
   
    %%% norm_lf
    cfg                 = [];
    cfg.sourcemodel.pos = pos135;
    cfg.grad            = fieldlinebeta2_head_fixedsensors;
    cfg.channel         = fieldlinebeta2_head_fixedsensors.label(1:128);
    cfg.headmodel       = headmodel_sphere;
    leadfield           = ft_prepare_leadfield(cfg);
    % 
    % mom         = mom135;
    % norm_mom    = norm(mom);    
    % orient      = mom ./ norm_mom;
    
    lf          = leadfield.leadfield{1,1} * orient;
    norm_lf  = vecnorm(lf);
  
    dmu3(l)     = sensor_noise_opm3(l) / norm_lf ;
      
end

%% dpu
figure;
plot(log10(sensor_noise_opm2), log10(dpu2), 'o');
hold on;
plot(log10(sensor_noise_opm), log10(dpu1), 'o');

xlabel('log10(noise)');
ylabel('log10(dpu)');


%%% fit a line to the first log-log plot.
coefficients1 = polyfit(log10(sensor_noise_opm2(3:6)), log10(dpu2(3:6)), 1); % Linear fit (degree 1)
a1 = coefficients1(1);
b1 = coefficients1(2);

tmp = linspace(1e-20, 1, 8);
fitted_line1 = a1 * log10(tmp) + b1;
plot(log10(tmp), fitted_line1, '--');

fitted_line_eq1 = sprintf('Monte-Carlo: y = %.2fx + %.2f', a1, b1);


%%% fit a line to the second log-log plot.
coefficients2 = polyfit(log10(sensor_noise_opm(3:6)), log10(dpu1(3:6)), 1); % Linear fit (degree 1)
a2 = coefficients2(1);
b2 = coefficients2(2);

fitted_line2 = a2 * log10(tmp) + b2;
plot(log10(tmp), fitted_line2, '--');

fitted_line_eq2 = sprintf('Jackknife: y = %.2fx + %.2f', a2, b2);




% Combine legends
legend('Monte-Carlo', 'Jackknife', fitted_line_eq1, fitted_line_eq2, 'Location', 'northwest');


%% dmu
figure;
plot(log10(sensor_noise_opm2), log10(dmu2), 'o');
hold on;
plot(log10(sensor_noise_opm), log10(dmu1), 'go');
hold on;
plot(log10(sensor_noise_opm3), log10(dmu3), 'o');

xlabel('log10(noise)');
ylabel('log10(dmu)');

figure;
plot(log10(sensor_noise_opm2), log10(dmu2_no_inflation), 'o');
hold on;
plot(log10(sensor_noise_opm), log10(dmu1_no_inflation), 'go');
hold on;
plot(log10(sensor_noise_opm3), log10(dmu3), 'o');


%% Publication figures

% fix units:
% 1. when we work in cm the leadfield units are in 1e2 T/(A*cm). So the dipole moment units should be 1e-2 A*cm so than the forward model gives the magnetic field in T. (see https://mailman.science.ru.nl/pipermail/fieldtrip/2014-September/008438.html)
% 2. go from A*cm -> A*m. I multiply by 1e-2
dmu2_no_inflation = 10^-4*dmu2_no_inflation;
dmu1_no_inflation = 10^-4*dmu1_no_inflation;
dmu3 = 10^-4* dmu3;

%%
figure
scatter(sensor_noise_opm2, dmu2_no_inflation, 150, 'o', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.7);
hold on
scatter(sensor_noise_opm2, dmu1_no_inflation, 150, 'd', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'g', 'MarkerFaceAlpha', 0.7);
hold on
scatter(sensor_noise_opm2, dmu3, 150, 's', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.7);

set(gca, 'XScale', 'log', 'YScale', 'log'); % Set the axes to log scale

grid on

% Fit lines in log-log space
% Fitting line for dmu2_no_inflation (blue)
p1 = polyfit(log10(sensor_noise_opm2), log10(dmu2_no_inflation), 1);
yfit1 = 10.^polyval(p1, log10(sensor_noise_opm2));
plot(sensor_noise_opm2, yfit1, '--b', 'LineWidth', 1.5);

% Fitting line for dmu1_no_inflation (green)
p2 = polyfit(log10(sensor_noise_opm2), log10(dmu1_no_inflation), 1);
yfit2 = 10.^polyval(p2, log10(sensor_noise_opm2));
plot(sensor_noise_opm2, yfit2, '--g', 'LineWidth', 1.5);

% Fitting line for dmu3 (red)
p3 = polyfit(log10(sensor_noise_opm2), log10(dmu3), 1);
yfit3 = 10.^polyval(p3, log10(sensor_noise_opm2));
plot(sensor_noise_opm2, yfit3, '--r', 'LineWidth', 1.5);

% Prepare equation strings for the legend in the form of log(y) = m*log(x) + c
eq1 = sprintf('log(y) = %.2f log(x) + %.2f', p1(1), p1(2));
eq2 = sprintf('log(y) = %.2f log(x) + %.2f', p2(1), p2(2));
eq3 = sprintf('log(y) = %.2f log(x) + %.2f', p3(1), p3(2));

% Add legend with fitted line equations
legend({'dmu2 (blue)', eq1, 'dmu1 (green)', eq2, 'dmu3 (red)', eq3}, 'Location', 'best');

xlabel('Sensor noise (T)', 'FontSize', 14);
ylabel('DMU_OPM (A*m)')

print('A2', '-dpng', '-r300'); % Save as PNG with 300 DPI resolution

