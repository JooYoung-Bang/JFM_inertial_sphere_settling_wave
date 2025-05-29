clc
close all
clear all

% example images are one of N1 settling results (Nylon 3/32", specific
% gravity of 1.12)

% read image
a=dir('*.tif'); 
num_img = numel(a); % # of tif images

for ii=1:num_img
    s = strcat('img (',num2str(ii),').tif'); 
    b = imread(s);
    img_3d(:,:,ii) = uint8(b(:,:)); %stack of images
end

% remove background
bg_img = min(img_3d,[],3); % background image
img_bgr = img_3d -bg_img;  % Background removed img (bgr:background removed)  

imshow(img_bgr(:,:,1))

%% Centroid Detection via Hough Transform (imfindcircles)
clc;
close all;

% Extract dimensions from the image stack
[num_x, num_z, num_img] = size(img_3d);

% Sampling parameters
freq = 20;            % Sampling frequency (Hz)
dt = 1 / freq;        % Sampling interval (s)

% Spatial resolution (pixels per millimeter)
px_mm = 714 / 100;    % Based on magnification

% Sphere radius in pixels (3/32 inch = 2.38 mm diameter)
radius_mm = 25.4 / 32 * 3 / 2;         % Radius in mm
radi_px = radius_mm * px_mm;          % Radius in pixels

% Radius range for circle detection (rounded to nearest integers)
r_min = radi_px * 1 - 2;
r_max = radi_px * 2;
r_range = [floor(r_min), floor(r_max)];

% Initialize tracking variables
num_loop = 0;         % Counter for frames with multiple detected candidates
centroid = nan(num_img, 2);   % Centroid coordinates [x, z]
prtc_dia = nan(num_img, 1);   % Particle diameter [px]

% Centroid estimation loop
for ii = 1:num_img
    frame = double(img_bgr(:, :, ii));  % Convert to double for processing
    
    % Detect circular features within the defined radius range
    [centers, radii] = imfindcircles(frame, r_range);
    
    if numel(radii) > 1
        % More than one candidate detected — select the largest
        num_loop = num_loop + 1;
        num_error(num_loop, 1) = ii;
        
        [~, idx] = max(radii);
        centroid(ii, :) = centers(idx, :);
        prtc_dia(ii, 1) = 2 * radii(idx);
    elseif numel(radii) == 1
        % Single candidate detected
        centroid(ii, :) = centers;
        prtc_dia(ii, 1) = 2 * radii;
    end
end
centroid(:,2) = -centroid(:,2);

%% Quality Check of Detection with Inset (Zoom-in View)
clc;
close all;

% Choose a frame to inspect detection quality
test_num = 1; 
rad_prtc = prtc_dia(test_num) / 2;

% Main figure: full frame
figure;
imshow(img_3d(:, :, test_num));
hold on;

% Overlay estimated centroid and fitted circle
scatter(centroid(test_num, 1), -centroid(test_num, 2), 'ro', 'filled');
viscircles([centroid(test_num, 1), -centroid(test_num, 2)], rad_prtc, 'LineWidth', 1);
title(['Detection Quality: Frame ', num2str(test_num)], 'Interpreter', 'latex');

% Define zoom window around centroid
zoom_size = 50; % Half-width of the zoom-in window in pixels
x_cent = round(centroid(test_num, 1));
y_cent = round(-centroid(test_num, 2));

% Ensure the window is within image bounds
x_min = max(x_cent - zoom_size, 1);
x_max = min(x_cent + zoom_size, size(img_3d, 2));
y_min = max(y_cent - zoom_size, 1);
y_max = min(y_cent + zoom_size, size(img_3d, 1));

% Create inset axes
inset_ax = axes('Position', [0.65, 0.65, 0.25, 0.25]); % [left bottom width height]
inset_img = img_3d(y_min:y_max, x_min:x_max, test_num);
imshow(inset_img, 'Parent', inset_ax);
hold(inset_ax, 'on');

% Overlay circle and centroid in the inset (adjust coordinates)
scatter(inset_ax, zoom_size + 1, zoom_size + 1, 'ro', 'filled');
viscircles(inset_ax, [zoom_size + 1, zoom_size + 1], rad_prtc, 'LineWidth', 1);

% Optional: Box on inset and axis off for cleaner display
axis(inset_ax, 'off');
box(inset_ax, 'on');

set(gcf,'position',[100 100 800 500])


%% Trajectory Filtering Using Savitzky-Golay Filter with data comparison (raw vs filtered)
clc
close all
num_filter = 5;     % Window size (must be odd)
ord_filter = 2;     % Polynomial order
filtered_centroid(:, 1) = sgolayfilt(centroid(:, 1), ord_filter, num_filter);
filtered_centroid(:, 2) = sgolayfilt(centroid(:, 2), ord_filter, num_filter);

% Velocity Estimation Using Finite Difference Scheme
% Note: Division by 10 converts velocity from mm/s to cm/s

velo_x_raw = (centroid(2:end, 1) - centroid(1:end-1, 1)) / px_mm / dt / 10;
velo_z_raw = (centroid(2:end, 2) - centroid(1:end-1, 2)) / px_mm / dt / 10;

velo_x_filtered = (filtered_centroid(2:end, 1) - filtered_centroid(1:end-1, 1)) / px_mm / dt / 10;
velo_z_filtered = (filtered_centroid(2:end, 2) - filtered_centroid(1:end-1, 2)) / px_mm / dt / 10;

% Time array centered between consecutive frames
tt_velo = (1:num_img-1) * dt + 0.5 * dt;


% Plot: Particle Trajectory (Raw and Filtered)
figure(1);
plot(centroid(:, 1), centroid(:, 2), '-bo');  % Raw trajectory
hold on;
plot(filtered_centroid(:, 1), filtered_centroid(:, 2), '-rx');  % Filtered trajectory
axis equal;
xlabel('$x^\prime$ (px)', 'Interpreter', 'latex');
ylabel('$z^\prime$ (px)', 'Interpreter', 'latex');
title('Particle Trajectory (Centroid Path)', 'Interpreter', 'latex');
legend('Raw', 'Filtered', 'Interpreter', 'latex');

figure(2);

subplot(2, 1, 1);
plot(tt_velo, velo_x_raw, '-bo');
hold on;
plot(tt_velo, velo_x_filtered, '-rx');
ylabel('$v^\prime_{x}$ (cm/s)', 'Interpreter', 'latex');
legend('Raw', 'Filtered', 'Interpreter', 'latex');
grid on;
title('Particle Velocity', 'Interpreter', 'latex');
subplot(2, 1, 2);
plot(tt_velo, velo_z_raw, '-bo');
hold on;
plot(tt_velo, velo_z_filtered, '-rx');
xlabel('$t$ (s)', 'Interpreter', 'latex');
ylabel('$v^\prime_{z}$ (cm/s)', 'Interpreter', 'latex');
grid on;

%% Load PIV Velocity Field and Define Wave & Particle Parameters

clc;
close all;

% Load pixel-to-physical unit conversion data
load flow_px_unit;

% Physical and fluid properties
prtc_dia_mm = 25.4 / 32 * 3;   % Particle diameter [mm]
rho_w = 1000;                  % Water density at 20°C [kg/m^3]
rho_p = 1120;                  % Particle density [kg/m^3]
Cm = 0.5;                      % Added mass coefficient for a sphere
nu = 1e-6;                     % Kinematic viscosity [m^2/s]
g = 9.81;                      % Gravitational acceleration [m/s^2]

% Density ratio
gamma = rho_p / rho_w;

% Define Wave Parameters (Case W1)

T = 0.657;              % Wave period [s]
h = 0.315;              % Water depth [m]
a = 0.0076;             % Wave amplitude [m]
omega = 2 * pi / T;     % Angular frequency [rad/s]

% Compute wavenumber and wavelength using linear wave theory
[k, lambda] = dispersion(h, omega);  % k: wavenumber [rad/m], lambda: wavelength [m]

% Wave-related dimensionless parameters
kh = k * h;                    % Relative depth
ka = k * a;                    % Wave steepness
eps = ka / tanh(kh);           % Nonlinearity parameter

% Synchronize Particle Centroid with PIV Timestamps
% Estimate particle position at PIV time steps via temporal midpoint
ctr_piv_matched = (filtered_centroid(1:end-1, :) + filtered_centroid(2:end, :)) / 2;

%% Estimate Undisturbed Flow Velocity at the Particle Position

clc;
close all;

% Dimensions of the PIV velocity field
[num_x, num_z, num_img] = size(uu_raw);
tt_dom = (1:num_img) * dt;  % Time domain [s]

% Grid spacing in the PIV domain (assumes uniform grid)
dx = xx_dom(2,1) - xx_dom(1,1);

% Threshold region of influence around the particle (in physical units)
% Removing PIV vectors within ~3 particle diameters
thres_affect = px_mm * prtc_dia_mm * 3;

% Initialize arrays for interpolated velocities
uu_est = nan(num_img, 1);
ww_est = nan(num_img, 1);

for kk = 1:num_img
    % Current frame index
    test_num = kk;
    
    % Initialize temporary fields
    xx_temp = xx_dom;
    zz_temp = zz_dom;
    uu_temp = uu_raw(:,:,test_num);
    ww_temp = ww_raw(:,:,test_num);
    
    % Mask velocities near the particle centroid
    for ii = 1:num_x
        for jj = 1:num_z
            dist_x = abs(xx_dom(ii,jj) - ctr_piv_matched(test_num,1));
            dist_z = abs(zz_dom(ii,jj) - ctr_piv_matched(test_num,2));
            if dist_x < thres_affect && dist_z < thres_affect 
                xx_temp(ii,jj) = NaN;
                zz_temp(ii,jj) = NaN;
                uu_temp(ii,jj) = NaN;
                ww_temp(ii,jj) = NaN;
            end
        end
    end

    % Restrict interpolation to the region below the particle release point
    zz_ini = 10; % Index corresponding to ~ 3a (3 x wave amplitude) below the wave trough
    xx_temp = double(xx_temp(:, zz_ini:end));
    zz_temp = double(zz_temp(:, zz_ini:end));
    uu_temp = double(uu_temp(:, zz_ini:end));
    ww_temp = double(ww_temp(:, zz_ini:end));

    % Logical mask of valid (non-NaN) data points
    valid_mask = ~isnan(xx_temp);
    
    
        
    % Interpolate horizontal velocity (u) at particle position
    uu_F = scatteredInterpolant(xx_temp(valid_mask), zz_temp(valid_mask), uu_temp(valid_mask), 'linear', 'linear');
    uu_est(kk,1) = uu_F([ctr_piv_matched(test_num,1), ctr_piv_matched(test_num,2)]);

    % Interpolate vertical velocity (w) at particle position
    ww_F = scatteredInterpolant(xx_temp(valid_mask), zz_temp(valid_mask), ww_temp(valid_mask), 'linear', 'linear');
    ww_est(kk,1) = ww_F([ctr_piv_matched(test_num,1), ctr_piv_matched(test_num,2)]);
    
end


%% Visualization of Results and Figures
clc;
close all;

% Convert velocities from pixels/frame to mm/s
ptr_velo = (filtered_centroid(2:end,:) - filtered_centroid(1:end-1,:)) / (dt * px_mm);  % Particle velocity (mm/s)
ptr_velo_fluid = [uu_est ww_est] / (dt * px_mm);  % Undisturbed fluid velocity (mm/s)

% Terminal settling velocity of the N1 particle (m/s)
terminal_velo = -0.0745;

% Define water surface position in pixels
z_surface = 102;

% Convert particle centroid coordinates from pixels to millimeters
ptr_centroid_mm(:,1) = ctr_piv_matched(:,1) / px_mm;
ptr_centroid_mm(:,2) = (z_surface + ctr_piv_matched(:,2)) / px_mm;

% Plot particle trajectory (in cm units)
figure(1);
scatter((ptr_centroid_mm(:,1) - ptr_centroid_mm(1,1)) / 10, ptr_centroid_mm(:,2) / 10, 'filled');
xlabel('$x$ (cm)', 'Interpreter', 'latex');
ylabel('$z$ (cm)', 'Interpreter', 'latex');
box on;
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
ylim([-27 0]);
set(gcf, 'Position', [100 100 400 800]);

% Plot particle and fluid velocities (in cm/s)
figure(2);

% Horizontal velocity
h_fig(1) = subplot(2,1,1);
plot(tt_dom, ptr_velo(:,1)/10, '-o');
hold on;
plot(tt_dom, ptr_velo_fluid(:,1)/10, '-o');
ylabel('$(\textrm{cm/s})$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

% Vertical velocity
h_fig(2) = subplot(2,1,2);
plot(tt_dom, ptr_velo(:,2)/10, '-o');
hold on;
plot(tt_dom, -ptr_velo_fluid(:,2)/10, '-o');
plot(tt_dom, -ptr_velo_fluid(:,2)/10 + terminal_velo * 100 * ones(num_img,1), '-o');
ylabel('$(\textrm{cm/s})$', 'Interpreter', 'latex');
xlabel('$t$ (s)', 'Interpreter', 'latex');
grid on
legend(h_fig(2), ...
    '$v_z$', '$u_z$', '$u_z - v_s$', ...
    'Interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 14);

legend(h_fig(1), ...
    '$v_x$', '$u_x$', ...
    'Interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 14);

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
set(gcf, 'Position', [500 100 600 600]);

% Plot slip velocities (in cm/s)
figure(3);

% Horizontal slip velocity
h_fig(1) = subplot(2,1,1);
plot(tt_dom, (ptr_velo(:,1) - ptr_velo_fluid(:,1)) / 10, '-*');
hold on;
plot(tt_dom, zeros(num_img,1), '-k');
plot(tt_dom, 0.05 * terminal_velo * 100 * ones(num_img,1), '--k');
plot(tt_dom, -0.05 * terminal_velo * 100 * ones(num_img,1), '--k');
ylabel('$(\textrm{cm/s})$', 'Interpreter', 'latex');
xlabel('$t$ (s)', 'Interpreter', 'latex');
legend('$v_x - u_x$', 'Interpreter', 'latex', 'Location', 'best');
grid on;

% Vertical slip velocity
h_fig(2) = subplot(2,1,2);
plot(tt_dom, (ptr_velo(:,2) + ptr_velo_fluid(:,2)) / 10, '-x');
hold on;
plot(tt_dom, terminal_velo * 100 * ones(num_img,1), '-k');
plot(tt_dom, 0.95 * terminal_velo * 100 * ones(num_img,1), '--k');
plot(tt_dom, 1.05 * terminal_velo * 100 * ones(num_img,1), '--k');
grid on
legend('$v_z - u_z$', '$-v_{s}$', 'Interpreter', 'latex', 'Location', 'best');
ylabel('$(\textrm{cm/s})$', 'Interpreter', 'latex');
xlabel('$t$ (s)', 'Interpreter', 'latex');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
set(gcf, 'Position', [1100 100 500 600]);

%% Dispersion function
function [k,lambda] = dispersion(h, omega)

% input h is water depth in metres
% input omega is angular frequency in in 1/s

% gravitational acceleration in m/s^2
g=9.81;

% solving dispersion relation omega^2 = g k tanh(kh)
om_nd = omega^2*h/g; % omega^2*(h/g)
f1 = @(x) x.*tanh(x) - om_nd; 
kh = fzero(f1,sqrt(om_nd));
k = kh/h; % wavenumber
lambda = 2*pi/k; % wavelength in m
end