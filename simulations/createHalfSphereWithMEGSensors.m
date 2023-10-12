function createHalfSphereWithMEGSensors(numSensors)
    % Create a half-sphere with MEG sensors on its edges
    
    % Define the radius of the sphere
    radius = 1;
    
    % Generate equally spaced points on the edges of a unit sphere
    theta = linspace(0, pi, numSensors);
    phi = linspace(0, 2 * pi, numSensors);
    
    % Create a grid of points on the half-sphere
    [Theta, Phi] = meshgrid(theta, phi);
    
    % Calculate the coordinates of the points on the half-sphere
    x = radius * sin(Theta) .* cos(Phi);
    y = radius * sin(Theta) .* sin(Phi);
    z = radius * cos(Theta);
    
    % Calculate sensor positions on the edges of the half-sphere
    sensorX = x(:);
    sensorY = y(:);
    sensorZ = z(:);
    
    % Plot the half-sphere
    figure;
    surf(x, y, z);
    axis equal;
    title(['Half-Sphere with ', num2str(numSensors), ' MEG sensors on edges']);
    
    % Plot MEG sensor positions as red points
    hold on;
    scatter3(sensorX, sensorY, sensorZ, 20, 'r', 'filled');
    hold off;
    
    % Label the MEG sensor points
    for i = 1:numel(sensorX)
        text(sensorX(i), sensorY(i), sensorZ(i), num2str(i), 'Color', 'k');
    end
end


