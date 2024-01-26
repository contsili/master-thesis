%% Single sphere head model
headmodel = [];
headmodel.r = 0.09;
headmodel.o = [0 0 0];
headmodel.unit = 'm';

[x,y,z] = sphere;
x = x * headmodel.r + headmodel.o(1);
y = y * headmodel.r + headmodel.o(2);
z = z * headmodel.r + headmodel.o(3);
h=surf(x,y,z);
set(h, 'FaceAlpha', 1);
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');

%% Sensor array

load fieldlinebeta2.mat; 
ft_plot_sens(fieldlinebeta2, 'label', 'no', 'axes', 0, 'orientation', 0) % looks good

%% Source model (points on a cortical surface)

% way 1
% cfg = [];
% cfg.method = 'basedoncortex';
% cfg.headshape = headmodel; % it needs a .fif file. What shall I do?
% sourcemodel = ft_prepare_sourcemodel(cfg)

% way 2: manually create points on my sphere
cd 'C:\Users\user\Documents\Courses\Internship\master-thesis\simulations\ni2-electrophys-master\private' % that's where mesh_sphere() is located

sourcemodel                       = mesh_sphere(5000, 'msphere');
sourcemodel                       = sourcemodel/11;
sourcemodel(sourcemodel(:,3)<0,:) = []; % keep the upper sphere only

% Plot the points on the surface of the sphere
figure;
h_sphere = surf(x, y, z);
hold on;
scatter3(sourcemodel(:,1), sourcemodel(:,2), sourcemodel(:,3), 2, 'r', 'filled'); % Adjust the size and color as needed
hold off;

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Points on the Surface of the Sphere');
legend('Sphere', 'Generated Points');

%% Leadfield

%% way 1: ft_prepare_leadfield() - I do not know what to put as "data" in the input which defines the dipole moments.
% cfg =[];
% cfg.sourcemodel.pos = sourcemodel;
% cfg.headmodel = headmodel;
% cfg.grad = fieldlinebeta2;
% data.grad = fieldlinebeta2;
% leadfield = ft_prepare_leadfield(cfg);

%% way 2: Define dipole moment manually
dipmom = [1 0 0]'; % I choose a random dipmom

lf144 = zeros(144, length(sourcemodel));
for i = 1: length(sourcemodel) 
    lf144(:,i) = ft_compute_leadfield(sourcemodel(i,:), fieldlinebeta2, headmodel) * dipmom; % you do not need to do it in a loop: lf144 = t_compute_leadfield(sourcemodel, fieldlinebeta2, headmodel) * dipmom, where sourcemodel: Ndip*3
end

norm_lf144 = vecnorm(lf144);
norm_lf144=  sqrt(norm_lf144); % Some high values and a lot low ones -> use logarithmic scale


%% way 3: I choose a tangential dipole always
lf144_tan = zeros(144, length(sourcemodel));

for i = 1:length(sourcemodel)

    angle =pi/2;

    % Convert to spherical coordinates
    [phi, theta, rho] = cart2sph(sourcemodel(i, 1), sourcemodel(i, 2), sourcemodel(i, 3));
    % Define the rotation axis as the phi vector
    rotation_axis = [sin(phi) -cos(phi) 0];
    % Define the rotation matrix using Rodrigues' rotation formula
    R = eye(3) + sin(angle) * cross_matrix(rotation_axis) + (1 - cos(angle)) * cross_matrix(rotation_axis)^2;
    % Rotate the point
    dipmom2(i,:) = (R * sourcemodel(i, :)')';
      
    lf144_tan(:, i) = ft_compute_leadfield(sourcemodel(i, :), fieldlinebeta2, headmodel) * dipmom2(i,:)';
end

norm_lf144_tan = vecnorm(lf144_tan);
norm_lf144_tan=  log(norm_lf144_tan); % Some high values and a lot low ones -> use logarithmic scale


% check if I picked the tangential dipoles
quiver3(sourcemodel(:, 1), sourcemodel(:, 2), sourcemodel(:, 3), ...
    dipmom2(:, 1), dipmom2(:, 2), dipmom2(:, 3), 'Color', 'b'); % looks good

% What if we turn the dipole to another tangential orientation (eg pointing always to the right (or left) and not to the top)?

%% way 4: I choose a normal dipole always
lf144_nor = zeros(144, length(sourcemodel));

for i = 1:length(sourcemodel)
    dipmom3(i,:) = sourcemodel(i, :);
    lf144_nor(:, i) = ft_compute_leadfield(sourcemodel(i, :), fieldlinebeta2, headmodel) * dipmom3(i,:)';
end

norm_lf144_nor = vecnorm(lf144_nor);
norm_lf144_nor=  log(norm_lf144_nor); % Some high values and a lot low ones -> use logarithmic scale

% check if I picked the tangential dipoles
quiver3(sourcemodel(:, 1), sourcemodel(:, 2), sourcemodel(:, 3), ...
    dipmom3(:, 1), dipmom3(:, 2), dipmom3(:, 3), 'Color', 'b'); % looks good

%% Plot the colorcoded sphere (color is ||L_opm||)- discretized colorcoding

figure;
scatter3(sourcemodel(:, 1), sourcemodel(:, 2), sourcemodel(:, 3), 50, norm_lf144, 'filled');
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Color-coded ||L_opm||');
grid on;

%% Plot the colorcoded sphere (color is ||L_opm||) - continuous colorcoding

figure;

% Triangulation for the sphere surface
trisurf(delaunay(sourcemodel(:, 1), sourcemodel(:, 2), sourcemodel(:, 3)), ...
    sourcemodel(:, 1), sourcemodel(:, 2), sourcemodel(:, 3), ...
    norm_lf144, 'FaceColor', 'interp', 'EdgeColor', 'none');

colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Color-coded norm-lf144-nor');
grid on;

hold on;
ft_plot_sens(fieldlinebeta2, 'label', 'no', 'axes', 0, 'orientation', 0) % , 'chantype', 'megmag') % looks good

%% Home-made functions
function C = cross_matrix(v) % used in the Rodrigues formula
    % This function returns the cross product matrix of a vector v
    C = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
end