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

% Without random sensor noise
leadfield1 = ni2_leadfield(sensors, headmodel, dippar1);
ni2_topoplot(sensors, leadfield1); colorbar % a reason why the topoplot is not good for opm_tangential_radial is because I have 2 signals from 1 location. But the lf thinks these are two different locations, so in some close locations the signals add up and some others they subtract leading to this highly variable topography. The solution would be to add up the signals that are picked up by the sensors in the same location: leadfield= leadfield2 + leadfield1 and then plot 261 locations, not 522. I do that in the next lines

% With random sensor noise
leadfield1 = ni2_leadfield(sensors, headmodel, dippar1, 0.1 * 10^-10);
ni2_topoplot(sensors, leadfield1); colorbar 


%% 2. add linearly opm_tangential and opm_radial

headmodel = ni2_headmodel('type', 'spherical', 'nshell', 1); 
sensors1 = ni2_sensors('type', 'opm_tangential');
sensors2 = ni2_sensors('type', 'opm_radial');

dippar1 = [0 0 6 1 0 0];

% With random sensor noise
leadfield1 = ni2_leadfield(sensors1, headmodel, dippar1);
leadfield2 = ni2_leadfield(sensors2, headmodel, dippar1);
leadfield= leadfield2 + leadfield1;

ni2_topoplot(sensors1, leadfield1); colorbar % opm_tangential
ni2_topoplot(sensors2, leadfield2); colorbar % opm_radial
ni2_topoplot(sensors1, leadfield); colorbar % opm_tangential + opm_radial


% Also I can add them up vectorially since the sensors have 90 deg angle
% difference:
leadfield= sqrt(leadfield2.^2 + leadfield1.^2);
ni2_topoplot(sensors1, leadfield); colorbar

%% 3. ctf275

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

%% TODO: Select different channels with cfg.channel, also find the typical sensor noise, also add correlatted noise (maybe this is added by adding dipoles inside the head (background brain activity) and outside the head (correlated noise from environmental sources like a car or elevator*))
% * a car or elevator moves regarding the MEG, so it is a moving source (we will model this later)

%% Select different number of channels 


% way1: in fieldtrip you can use ft_channelselection and the 'gui' option. But you have to do this manually one by one channel


% way2: clustering algorithm to reduce the number of sensors while keeping
% them homogeneously distributed. But it does not inherently ensure that the new sensor positions (cluster centroids) will lie on the sphere
sensors1 = ni2_sensors('type', 'opm_radial');

% Assuming sensors1.coilpos is your sensor positions matrix
sensor_positions = sensors1.coilpos;

% Define the number of clusters (i.e., the reduced number of sensors you want)
num_clusters = 29; % change this to your desired number

% Run k-means clustering
[cluster_idx, cluster_centroid] = kmeans(sensor_positions, num_clusters);

% Now, cluster_centroid contains the positions of the reduced sensors
% sensors1.coilpos = cluster_centroid;
 scatter3(cluster_centroid(:,1),cluster_centroid(:,2),cluster_centroid(:,3), 'r')




% way3:  computing the Euclidean distance between each pair of sensors. It then enters a loop where in each iteration it finds and removes the sensor that has the closest neighbor. 
sensors1 = ni2_sensors('type', 'opm_radial');

% Assuming sensors1.coilpos is your sensor positions matrix
sensor_positions = sensors1.coilpos;

% Assuming 'sensor_positions' is your 261x3 sensor positions matrix
num_sensors = size(sensor_positions, 1);

% Define the desired number of sensors
desired_num_sensors = 29; % change this to your desired number

% Compute the Euclidean distance between each pair of sensors
distances = pdist2(sensor_positions, sensor_positions);

% Set the diagonal to infinity so that a sensor is not considered its own neighbor
distances(1:num_sensors+1:end) = Inf;

while num_sensors > desired_num_sensors
    % Find the sensor that has the closest neighbor
    [~, sensor_to_remove] = min(min(distances));
    
    % Remove this sensor from the distance matrix
    distances(sensor_to_remove, :) = [];
    distances(:, sensor_to_remove) = [];
    
    % Remove this sensor from 'sensor_positions'
    sensor_positions(sensor_to_remove, :) = [];
    
    % Update 'num_sensors'
    num_sensors = num_sensors - 1;
end

% Now, 'sensor_positions' contains the positions of the reduced sensors
reduced_sensor_positions1 = sensor_positions;
figure; scatter3(reduced_sensor_positions1(:,1),reduced_sensor_positions1(:,2),reduced_sensor_positions1(:,3))




% way4: Spherical distance

sensors1 = ni2_sensors('type', 'opm_radial');

% Assuming sensors1.coilpos is your sensor positions matrix
sensor_positions = sensors1.coilpos;

% Convert Cartesian coordinates to spherical coordinates
[azimuth,elevation,~] = cart2sph(sensor_positions(:,1),sensor_positions(:,2),sensor_positions(:,3));

% Convert radians to degrees
lat = rad2deg(elevation);
lon = rad2deg(azimuth);

% Combine latitude and longitude into one matrix
sensor_positions_spherical = [lat lon];


% Assuming 'sensor_positions' is your 261x3 sensor positions matrix
num_sensors = size(sensor_positions_spherical, 1);

% Define the desired number of sensors
desired_num_sensors = 29; % change this to your desired number

% Compute the spherical distance between each pair of sensors
distances = zeros(num_sensors);
for i = 1:num_sensors
    for j = i+1:num_sensors
        lat1 = sensor_positions_spherical(i,1);
        lon1 = sensor_positions_spherical(i,2);
        lat2 = sensor_positions_spherical(j,1);
        lon2 = sensor_positions_spherical(j,2);
        
        R = 1; 
        dlat = deg2rad(lat2-lat1);
        dlon = deg2rad(lon2-lon1); 
        a = sin(dlat/2) * sin(dlat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2) * sin(dlon/2); % Haversine formula
        c = 2 * atan2(sqrt(a), sqrt(1-a));    

        d = R * c; 
        
        distances(i,j) = d;
        distances(j,i) = d;
    end
end

% Set the diagonal to infinity so that a sensor is not considered its own neighbor
distances(1:num_sensors+1:end) = Inf;

while num_sensors > desired_num_sensors
    % Find the sensor that has the closest neighbor
    [~, sensor_to_remove] = min(min(distances));
    
    % Remove this sensor from the distance matrix
    distances(sensor_to_remove, :) = [];
    distances(:, sensor_to_remove) = [];
    
    % Remove this sensor from 'sensor_positions'
    sensor_positions_spherical(sensor_to_remove, :) = [];
    
    % Update 'num_sensors'
    num_sensors = num_sensors - 1;
end

% Now, 'sensor_positions' contains the positions of the reduced sensors
reduced_sensor_positions = sensor_positions_spherical;
% Convert spherical coordinates back to Cartesian for plotting
[x, y, z] = sph2cart(deg2rad(reduced_sensor_positions(:,2)), deg2rad(reduced_sensor_positions(:,1)), 1);

% Now you can use scatter3 to plot the points
figure; scatter3(x, y, z);






%way5: it seems like many sensors are gathered in the low parts of the head
%(z close to zero). a more homogeneous distribution of sensors around the head, you could modify the sensor reduction process to take into account the z-coordinate (height) of the sensors. One way to do this is to divide the sensors into different height levels and ensure that each level has a minimum number of sensors. 

sensors1 = ni2_sensors('type', 'opm_radial');

% Assuming sensors1.coilpos is your sensor positions matrix
sensor_positions = sensors1.coilpos;
% Assuming 'sensor_positions' is your 261x3 sensor positions matrix
num_sensors = size(sensor_positions, 1);

% Define the desired number of sensors
desired_num_sensors = 32; % change this to your desired number

% Compute the spherical distance between each pair of sensors
distances = zeros(num_sensors);
for i = 1:num_sensors
    for j = i+1:num_sensors
        lat1 = sensor_positions(i,1);
        lon1 = sensor_positions(i,2);
        lat2 = sensor_positions(j,1);
        lon2 = sensor_positions(j,2);
        
        R = 6371; % Radius of the earth in km
        dlat = deg2rad(lat2-lat1);
        dlon = deg2rad(lon2-lon1); 
        a = sin(dlat/2) * sin(dlat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2) * sin(dlon/2); 
        c = 2 * atan2(sqrt(a), sqrt(1-a)); 
        d = R * c; % Distance in km
        
        distances(i,j) = d;
        distances(j,i) = d;
    end
end

% Set the diagonal to infinity so that a sensor is not considered its own neighbor
distances(1:num_sensors+1:end) = Inf;

% Divide sensors into different height levels based on z-coordinate
num_levels = 10; % change this to your desired number of levels
[~,~,level_indices] = histcounts(sensor_positions(:,3), num_levels);

while num_sensors > desired_num_sensors
    % Find the level with the most sensors
    level_counts = histcounts(level_indices, num_levels);
    [~, level_with_most_sensors] = max(level_counts);
    
    % Find the sensors in this level
    sensors_in_level = find(level_indices == level_with_most_sensors);
    
    % Find the sensor in this level that has the closest neighbor
    [~, sensor_to_remove_in_level] = min(min(distances(sensors_in_level, sensors_in_level)));
    sensor_to_remove = sensors_in_level(sensor_to_remove_in_level);
    
    % Remove this sensor from the distance matrix
    distances(sensor_to_remove, :) = [];
    distances(:, sensor_to_remove) = [];
    
    % Remove this sensor from 'sensor_positions'
    sensor_positions(sensor_to_remove, :) = [];
    
    % Update 'num_sensors' and 'level_indices'
    num_sensors = num_sensors - 1;
    level_indices(sensor_to_remove) = [];
end

% Now, 'sensor_positions' contains the positions of the reduced sensors
reduced_sensor_positions2 = sensor_positions;

figure; scatter3(reduced_sensor_positions2(:,1),reduced_sensor_positions2(:,2),reduced_sensor_positions2(:,3))





%way6: make less refinements in cosahedron
sensors1 = ni2_sensors('type', 'opm_tangential');
figure; ft_plot_sens(sensors1, 'elecshape', 'disc', 'elecsize', 1, 'label', 'off'); % 46 sensors




% way7:longitude lines
sensors1 = ni2_sensors('type', 'opm_radial');

% Assuming sensors1.coilpos is your sensor positions matrix
sensor_positions = sensors1.coilpos;

% Assuming 'sensor_positions' is your 261x3 sensor positions matrix
num_sensors = size(sensor_positions, 1);

% Define the desired number of sensors
desired_num_sensors = 32; % change this to your desired number

% Compute the spherical distance between each pair of sensors
distances = zeros(num_sensors);
for i = 1:num_sensors
    for j = i+1:num_sensors
        lat1 = sensor_positions(i,1);
        lon1 = sensor_positions(i,2);
        lat2 = sensor_positions(j,1);
        lon2 = sensor_positions(j,2);
        
        R = 6371; % Radius of the earth in km
        dlat = deg2rad(lat2-lat1);
        dlon = deg2rad(lon2-lon1); 
        a = sin(dlat/2) * sin(dlat/2) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2) * sin(dlon/2); 
        c = 2 * atan2(sqrt(a), sqrt(1-a)); 
        d = R * c; % Distance in km
        
        distances(i,j) = d;
        distances(j,i) = d;
    end
end

% Set the diagonal to infinity so that a sensor is not considered its own neighbor
distances(1:num_sensors+1:end) = Inf;

% Divide sensors into different longitude isolines based on x-coordinate
num_isolines = 15; % change this to your desired number of isolines
[~,~,isolines_indices] = histcounts(sensor_positions(:,1), num_isolines);

while num_sensors > desired_num_sensors
    % Find the isoline with the most sensors
    isoline_counts = histcounts(isolines_indices, num_isolines);
    [~, isoline_with_most_sensors] = max(isoline_counts);
    
    % Find the sensors in this isoline
    sensors_in_isoline = find(isolines_indices == isoline_with_most_sensors);
    
    % Find the sensor in this isoline that has the closest neighbor
    [~, sensor_to_remove_in_isoline] = min(min(distances(sensors_in_isoline, sensors_in_isoline)));
    sensor_to_remove = sensors_in_isoline(sensor_to_remove_in_isoline);
    
    % Remove this sensor from the distance matrix
    distances(sensor_to_remove, :) = [];
    distances(:, sensor_to_remove) = [];
    
    % Remove this sensor from 'sensor_positions'
    sensor_positions(sensor_to_remove, :) = [];
    
    % Update 'num_sensors' and 'isolines_indices'
    num_sensors = num_sensors - 1;
    isolines_indices(sensor_to_remove) = [];
end

% Now, 'sensor_positions' contains the positions of the reduced sensors
reduced_sensor_positions3 = sensor_positions;
figure; scatter3(reduced_sensor_positions3(:,1),reduced_sensor_positions3(:,2),reduced_sensor_positions3(:,3))

