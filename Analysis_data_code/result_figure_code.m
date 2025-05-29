clc;
close all;
clear all;

% this code includes how I draw figure 5 - 9. 

% This script plots the normalized particle trajectories in the x-z plane
% from experimental data for different types of particles.
% Normalized velocities (v'_s/c_w) for wave case W1 is provided (e.g.,
% N1_vs_norm, N3_vs_norm,...)

% Each particle dataset (e.g., nylon_332, nylon_18, etc.) contains 10 experiments (cells).
% Column structure:
% 1st: time [s]
% 2nd-3rd: position x, z [mm]
% 4th-5th: particle velocity v_x, v_z [mm/s]
% 6th-7th: undisturbed flow velocity u_x, u_z at particle position [mm/s]

load total_wave_filt  % Load processed particle trajectory data
%% Parameters (W1 wave case)
h_dep = 315;           % Water depth in mm
wave_length = 670.6;   % Wavelength in mm

%% Figure 5: Normalized Particle Trajectories in x-z Plane
% The initial x-position of each trajectory is shifted to zero
% for easy visual comparison. 

clc
close all
figure(1);

% Subplot (a): Nylon 332
h(1) = subplot(1, 5, 1);
hold on;
for kk = 1:10
    x_norm = (nylon_332{1, kk}(:, 2) - nylon_332{1, kk}(1, 2)) / wave_length;
    z_norm = nylon_332{1, kk}(:, 3) / h_dep;
    plot(x_norm, z_norm, '-x', 'MarkerSize', 6);
end
box on;
xlabel('$x/L$', 'Interpreter', 'latex');
ylabel('$z/h$', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
xlim([-0.04 0.04]);
ylim([-0.85 0]);

% Subplot (b): Nylon 18
h(2) = subplot(1, 5, 2);
hold on;
for kk = 1:10
    x_norm = (nylon_18{1, kk}(:, 2) - nylon_18{1, kk}(1, 2)) / wave_length;
    z_norm = nylon_18{1, kk}(:, 3) / h_dep;
    plot(x_norm, z_norm, '-x', 'MarkerSize', 6);
end
box on;
xlabel('$x/L$', 'Interpreter', 'latex');
grid on;
xlim([-0.04 0.04]);
ylim([-0.85 0]);

% Subplot (c): Nylon 14
h(3) = subplot(1, 5, 3);
hold on;
for kk = 1:10
    x_norm = (nylon_14{1, kk}(:, 2) - nylon_14{1, kk}(1, 2)) / wave_length;
    z_norm = nylon_14{1, kk}(:, 3) / h_dep;
    plot(x_norm, z_norm, '-x', 'MarkerSize', 6);
end
box on;
xlabel('$x/L$', 'Interpreter', 'latex');
grid on;
xlim([-0.07 0.07]);
ylim([-0.85 0]);
xticks([-0.04 0 0.04]);

% Subplot (d): Nylon 12 6mm
h(4) = subplot(1, 5, 4);
hold on;
for kk = 1:10
    x_norm = (nylon12_6{1, kk}(:, 2) - nylon12_6{1, kk}(1, 2)) / wave_length;
    z_norm = nylon12_6{1, kk}(:, 3) / h_dep;
    plot(x_norm, z_norm, '-x', 'MarkerSize', 6);
end
box on;
xlabel('$x/L$', 'Interpreter', 'latex');
grid on;
xlim([-0.04 0.04]);
ylim([-0.85 0]);

% Subplot (e): Torlon 18
h(5) = subplot(1, 5, 5);
hold on;
for kk = 1:10
    x_norm = (torlon_18{1, kk}(:, 2) - torlon_18{1, kk}(1, 2)) / wave_length;
    z_norm = torlon_18{1, kk}(:, 3) / h_dep;
    plot(x_norm, z_norm, '-x', 'MarkerSize', 6);
end
box on;
xlabel('$x/L$', 'Interpreter', 'latex');
grid on;
xlim([-0.04 0.04]);
ylim([-0.85 0]);

% Adjust subplot layout
set(gcf, 'Position', [100, 100, 1300, 550]);  % Set figure window size

sz_x = 0.125;      % Width of each subplot
intv_x = 0.184;    % Horizontal interval between subplots

% Manually set subplot positions for consistent spacing
for i = 1:5
    set(h(i), 'Position', [0.08 + (i-1)*intv_x, 0.13, sz_x, 0.8]);
end

% Add subplot labels (a–e)
text(h(1), -0.037, -0.03, '(a)', 'FontSize', 17, 'Interpreter', 'latex');
text(h(2), -0.037, -0.03, '(b)', 'FontSize', 17, 'Interpreter', 'latex');
text(h(3), -0.065, -0.03, '(c)', 'FontSize', 17, 'Interpreter', 'latex');
text(h(4), -0.037, -0.03, '(d)', 'FontSize', 17, 'Interpreter', 'latex');
text(h(5), -0.037, -0.03, '(e)', 'FontSize', 17, 'Interpreter', 'latex');

% Set global font size for all elements
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);

% Export options (uncomment to save)
% saveas(gcf, 'entire_trajectory_normalized.jpg');
% saveas(gcf, 'entire_trajectory_normalized.fig');
% saveas(gcf, 'entire_trajectory_normalized.eps', 'epsc');


%% Figure 6: Representative Velocity Profiles

clc
close all

% Select specific experiment numbers for Nylon and Torlon particles
test_N1 = 1;
test_T1 = 1;

figure(1);

% Subplot (a): Horizontal velocity components (v_x and u_x) for Nylon 332
h(1) = subplot(2, 2, 1);
plot(nylon_332{1, test_N1}(:, 1) * omega, nylon_332{1, test_N1}(:, 4) / c_w, '-s');
hold on;
plot(nylon_332{1, test_N1}(:, 1) * omega, nylon_332{1, test_N1}(:, 6) / c_w, '-x');
ylim([-0.06 0.06]);
yticks(-0.06:0.03:0.06);
xlim([0 35]);
grid on;
legend('$v_x$', '$u_x$', 'Interpreter', 'latex');

% Subplot (c): Vertical velocity components (v_z, u_z, u_z - v_s) for Nylon 332
h(3) = subplot(2, 2, 3);
plot(nylon_332{1, test_N1}(:, 1) * omega, nylon_332{1, test_N1}(:, 5) / c_w, '-s');
hold on;
plot(nylon_332{1, test_N1}(:, 1) * omega, nylon_332{1, test_N1}(:, 7) / c_w, '-o');
plot(nylon_332{1, test_N1}(:, 1) * omega, nylon_332{1, test_N1}(:, 7) / c_w + N1_vs_norm, '-x');
xlabel('$t$', 'Interpreter', 'latex');
ylim([-0.15 0.05]);
yticks(-0.15:0.05:0.05);
xlim([0 35]);
grid on;
legend('$v_z$', '$u_z$', '$u_z - v_s$', 'Interpreter', 'latex', 'NumColumns', 3, 'Location', 'south');

% Subplot (b): Horizontal velocity components (v_x and u_x) for Torlon 18
h(2) = subplot(2, 2, 2);
plot(torlon_18{1, test_T1}(:, 1) * omega, torlon_18{1, test_T1}(:, 4) / c_w, '-s');
hold on;
plot(torlon_18{1, test_T1}(:, 1) * omega, torlon_18{1, test_T1}(:, 6) / c_w, '-x');
ylim([-0.06 0.06]);
yticks(-0.06:0.03:0.06);
grid on;
legend('$v_x$', '$u_x$', 'Interpreter', 'latex');

% Subplot (d): Vertical velocity components (v_z, u_z, u_z - v_s) for Torlon 18
h(4) = subplot(2, 2, 4);
plot(torlon_18{1, test_T1}(:, 1) * omega, torlon_18{1, test_T1}(:, 5) / c_w, '-s');
hold on;
plot(torlon_18{1, test_T1}(:, 1) * omega, torlon_18{1, test_T1}(:, 7) / c_w, '-o');
plot(torlon_18{1, test_T1}(:, 1) * omega, torlon_18{1, test_T1}(:, 7) / c_w + T1_vs_norm, '-x');
xlabel('$t$', 'Interpreter', 'latex');
ylim([-0.3 0.05]);
yticks(-0.3:0.05:0.05);
grid on;
legend('$v_z$', '$u_z$', '$u_z - v_s$', 'Interpreter', 'latex', 'NumColumns', 3, 'Location', 'south');

% Adjust font sizes and figure dimensions
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
set(gcf, 'Position', [100, 100, 800, 550]);

% Add subplot labels (a–d)
text(h(1), -8.5, 0.06, '(a)', 'FontSize', 14, 'Interpreter', 'latex');
text(h(2), -4, 0.06, '(b)', 'FontSize', 14, 'Interpreter', 'latex');
text(h(3), -8.5, 0.05, '(c)', 'FontSize', 14, 'Interpreter', 'latex');
text(h(4), -4, 0.05, '(d)', 'FontSize', 14, 'Interpreter', 'latex');

% Export options (uncomment to save)
% saveas(gcf, 'Representative_velo_N1_T1.jpg');
% saveas(gcf, 'Representative_velo_N1_T1.fig');
% saveas(gcf, 'Representative_velo_N1_T1.eps', 'epsc');

%% Figure 7: Relative velocities in waves
clc
close all

lin_wid1 = 0.5;
lin_wid2 = 1.3;

figure(1)

% --- Nylon 33.2 (N1)
h(1) = subplot(5,2,1);
hold on
for kk = 1:10
    plot(nylon_332{1,kk}(:,1)*omega, ...
        (nylon_332{1,kk}(:,4) - nylon_332{1,kk}(:,6)) / c_w, '-x');
end
box on
plot([0 35], 0.05*N1_vs_norm*ones(2,1), '--k', 'linewidth', lin_wid1);
plot([0 35], -0.05*N1_vs_norm*ones(2,1), '--k', 'linewidth', lin_wid1);
plot([0 35], [0 0], '-k', 'linewidth', lin_wid2);
xlim([0 35]);
ylim([-0.06 0.06]);
grid on
ylabel('$(v-u)_x$', 'Interpreter', 'latex');

h(2) = subplot(5,2,2);
hold on
for kk = 1:10
    plot(nylon_332{1,kk}(:,1)*omega, ...
        (nylon_332{1,kk}(:,5) - nylon_332{1,kk}(:,7)) / c_w, '-x');
end
plot(omega * tau_pSN(1) * 4 * ones(2,1), [-0.08 0.05], '--r', 'linewidth', lin_wid2);
box on
plot([0 35], 1.05*N1_vs_norm*ones(2,1), '--k', 'linewidth', lin_wid1);
plot([0 35], 0.95*N1_vs_norm*ones(2,1), '--k', 'linewidth', lin_wid1);
plot([0 35], N1_vs_norm*ones(2,1), '-k', 'linewidth', lin_wid2);
xlim([0 35]);
ylim([-0.08 0.05]);
yticks([-0.05 0 0.05]);
grid on
ylabel('$(v-u)_z$', 'Interpreter', 'latex');

% --- Nylon 18 (N2)
h(3) = subplot(5,2,3);
hold on
for kk = 1:10
    plot(nylon_18{1,kk}(:,1)*omega, ...
        (nylon_18{1,kk}(:,4) - nylon_18{1,kk}(:,6)) / c_w, '-x');
end
box on
plot([0 30], 0.05*N2_vs_norm*ones(2,1), '--k', 'linewidth', lin_wid1);
plot([0 30], -0.05*N2_vs_norm*ones(2,1), '--k', 'linewidth', lin_wid1);
plot([0 30], [0 0], '-k', 'linewidth', lin_wid2);
grid on
ylim([-0.06 0.06]);
ylabel('$(v-u)_x$', 'Interpreter', 'latex');

h(4) = subplot(5,2,4);
hold on
for kk = 1:10
    plot(nylon_18{1,kk}(:,1)*omega, ...
        (nylon_18{1,kk}(:,5) - nylon_18{1,kk}(:,7)) / c_w, '-x');
end
box on
plot(omega * tau_pSN(2) * 4 * ones(2,1), [-0.1 0.05], '--r', 'linewidth', lin_wid2);
plot([0 30], [0.95 1.05]*N2_vs_norm, '--k', 'linewidth', lin_wid1);
plot([0 30], N2_vs_norm*ones(2,1), '-k', 'linewidth', lin_wid2);
grid on
ylim([-0.1 0.05]);
ylabel('$(v-u)_z$', 'Interpreter', 'latex');

% --- Nylon 14 (N3)
h(5) = subplot(5,2,5);
hold on
for kk = 1:10
    plot(nylon_14{1,kk}(:,1)*omega, ...
        (nylon_14{1,kk}(:,4) - nylon_14{1,kk}(:,6)) / c_w, '-x');
end
box on
plot([0 20], [0.05 -0.05]*N3_vs_norm, '--k', 'linewidth', lin_wid1);
plot([0 20], [0 0], '-k', 'linewidth', lin_wid2);
grid on
ylim([-0.06 0.06]);
ylabel('$(v-u)_x$', 'Interpreter', 'latex');

h(6) = subplot(5,2,6);
hold on
for kk = 1:10
    plot(nylon_14{1,kk}(:,1)*omega, ...
        (nylon_14{1,kk}(:,5) - nylon_14{1,kk}(:,7)) / c_w, '-x');
end
box on
plot([0 20], [0.95 1.05]*N3_vs_norm, '--k', 'linewidth', lin_wid1);
plot([0 20], N3_vs_norm*ones(2,1), '-k', 'linewidth', lin_wid2);
plot(omega * tau_pSN(3) * 4 * ones(2,1), [-0.16 0.05], '--r', 'linewidth', lin_wid2);
grid on
ylim([-0.16 0.05]);
yticks([-0.15 -0.1 -0.05 0 0.05]);
ylabel('$(v-u)_z$', 'Interpreter', 'latex');

% --- Nylon 12.6 (N4)
h(7) = subplot(5,2,7);
hold on
for kk = 1:10
    plot(nylon12_6{1,kk}(:,1)*omega, ...
        (nylon12_6{1,kk}(:,4) - nylon12_6{1,kk}(:,6)) / c_w, '-x');
end
box on
plot([0 45], [0.05 -0.05]*N4_vs_norm, '--k', 'linewidth', lin_wid1);
plot([0 45], [0 0], '-k', 'linewidth', lin_wid2);
xlim([0 45]);
grid on
ylim([-0.06 0.06]);
ylabel('$(v-u)_x$', 'Interpreter', 'latex');

h(8) = subplot(5,2,8);
hold on
for kk = 1:10
    plot(nylon12_6{1,kk}(:,1)*omega, ...
        (nylon12_6{1,kk}(:,5) - nylon12_6{1,kk}(:,7)) / c_w, '-x');
end
box on
plot([0 45], [0.95 1.05]*N4_vs_norm, '--k', 'linewidth', lin_wid1);
plot([0 45], N4_vs_norm*ones(2,1), '-k', 'linewidth', lin_wid2);
plot(omega * tau_pSN(4) * 4 * ones(2,1), [-0.07 0.05], '--r', 'linewidth', lin_wid2);
xlim([0 45]);
grid on
ylim([-0.07 0.05]);
yticks([-0.05 0 0.05]);
ylabel('$(v-u)_z$', 'Interpreter', 'latex');

% --- Torlon 18 (T1)
h(9) = subplot(5,2,9);
hold on
for kk = 1:10
    plot(torlon_18{1,kk}(:,1)*omega, ...
        (torlon_18{1,kk}(:,4) - torlon_18{1,kk}(:,6)) / c_w, '-x');
end
box on
plot([0 15], [0.05 -0.05]*T1_vs_norm, '--k', 'linewidth', lin_wid1);
plot([0 15], [0 0], '-k', 'linewidth', lin_wid2);
ylim([-0.06 0.06]);
ylabel('$(v-u)_x$', 'Interpreter', 'latex');
xlabel('$t$', 'Interpreter', 'latex');
grid on

h(10) = subplot(5,2,10);
hold on
for kk = 1:10
    plot(torlon_18{1,kk}(:,1)*omega, ...
        (torlon_18{1,kk}(:,5) - torlon_18{1,kk}(:,7)) / c_w, '-x');
end
box on
plot([0 15], [0.95 1.05]*T1_vs_norm, '--k', 'linewidth', lin_wid1);
plot([0 15], T1_vs_norm*ones(2,1), '-k', 'linewidth', lin_wid2);
plot(omega * tau_pSN(5) * 4 * ones(2,1), [-0.185 0.05], '--r', 'linewidth', lin_wid2);
ylim([-0.185 0.05]);
yticks([-0.15 -0.1 -0.05 0 0.05]);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$(v-u)_z$', 'Interpreter', 'latex');
grid on

% --- Subplot Labels
text_ft_sz = 15;
x_pos = -1.4;
y_pos1 = 0.06;
y_pos2 = 0.05;

labels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)'};
x_offsets = [7 7 6 6 4 4 9 9 3 3];
y_positions = [y_pos1 y_pos2 y_pos1 y_pos2 y_pos1 y_pos2 y_pos1 y_pos2 y_pos1 y_pos2];

for i = 1:10
    text(h(i), x_pos*x_offsets(i), y_positions(i), labels{i}, ...
        'FontSize', text_ft_sz, 'FontWeight', 'bold', 'Interpreter', 'latex');
end

set(findall(gcf,'-property','MarkerSize'),'MarkerSize',5)
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% subplot positioning
set(gcf,'position',[100 0 800 1600])
set(h(1),'position',[0.13 0.82 0.35 0.14])
set(h(3),'position',[0.13 0.63 0.35 0.14])
set(h(5),'position',[0.13 0.44 0.35 0.14])
set(h(7),'position',[0.13 0.25 0.35 0.14])
set(h(9),'position',[0.13 0.06 0.35 0.14])

set(h(2),'position',[0.60 0.82 0.35 0.14])
set(h(4),'position',[0.60 0.63 0.35 0.14])
set(h(6),'position',[0.60 0.44 0.35 0.14])
set(h(8),'position',[0.60 0.25 0.35 0.14])
set(h(10),'position',[0.60 0.06 0.35 0.14])

% % Optional: Save Figure
% saveas(figure(1), 'Entire_Relative_velocities.jpg');
% saveas(figure(1), 'Entire_Relative_velocities.fig');
% saveas(figure(1), 'Entire_Relative_velocities.eps', 'epsc');

%% Figure 8: Evolution of (v - u + v_s) for N1 particle

clc
close all

% Particle parameters
dp = 25.4/32 * 3;            % N1 particle diameter [mm] (3/32 inch)
vs_N1_dim = 74.5;            % terminal settling velocity v'_s [mm/s]
vs_N1 = 0.073;               % dimensionless terminal settling velocity v_s/c_w
omega = 9.56;                % Angular frequency
St_N1 = omega * tau_pSN(1);  % Stokes number for N1
Re_pt = vs_N1_dim * dp;      % Particle Reynolds number estimate

gamma = 1.12;                % Empirical coefficient

% Identify onset index tt_ini based on vertical slip exceeding threshold
for kk = 1:10
    for jj = 1:length(nylon_332{kk})
        slip_x = nylon_332{kk}(jj,4) - nylon_332{kk}(jj,6);
        slip_z = nylon_332{kk}(jj,5) - nylon_332{kk}(jj,7);
        if abs(slip_z) > 0.7 * vs_N1_dim
            tt_ini(kk,1) = jj;
            break
        end
    end
end

% Calculate instantaneous Re_p and St based on initial slip
for kk = 1:10
    jj = tt_ini(kk);
    slip_x = nylon_332{kk}(jj,4) - nylon_332{kk}(jj,6);
    slip_z = nylon_332{kk}(jj,5) - nylon_332{kk}(jj,7);
    Re_p(kk) = sqrt(slip_x^2 + slip_z^2) * dp;
    St(kk) = dp^2 * (gamma + 0.5) / (18 * (1 + 0.15 * Re_p(kk)^0.687)) * omega;
end

% Plotting parameters
lin_wid = 1.5;
x_max = 20;
x_min = 0;

% Manually identified t_x0 indices (based on slip_x behavior)
tt_ini_x = [5; 6; 6; 5; 7; 7; 5; 7; 7; 5];

clc
close all
figure(1)

% Subplot (a): Normalized horizontal slip decay
h(1) = subplot(1,2,1);
hold on

for kk = 1:10
    tt = nylon_332{kk}(:,1) * omega;
    slip_u = (nylon_332{kk}(:,4) - nylon_332{kk}(:,6)) / c_w;
    t0 = tt_ini_x(kk);
    plot(tt(t0:end) - tt(t0), slip_u(t0:end) / slip_u(t0), '-x');
end

tt_exp = 0:0.1:30;
plot(tt_exp, exp(-tt_exp / St_N1), '-k', 'LineWidth', lin_wid);
box on
xlabel('$(t - t_{x0})$', 'Interpreter', 'latex');
ylabel('$(v - u)_x / (v - u)_{x,t_{x0}}$', 'Interpreter', 'latex');
ylim([-1 1])
xlim([x_min x_max])

% Subplot (b): Normalized vertical slip + vs decay
h(2) = subplot(1,2,2);
hold on

for kk = 1:10
    tt = nylon_332{kk}(:,1) * omega;
    slip_w = (nylon_332{kk}(:,5) - nylon_332{kk}(:,7)) / c_w;
    t0 = tt_ini(kk);
    plot(tt(t0:end) - tt(t0), (slip_w(t0:end) + vs_N1) / (slip_w(t0) + vs_N1), '-x');
end

plot(tt_exp, exp(-tt_exp / St_N1), '-k', 'LineWidth', lin_wid);
box on
xlabel('$(t - t_{z0})$', 'Interpreter', 'latex');
ylabel('$((v - u)_z + v_s) / ((v - u)_{z,t_{z0}} + v_s)$', 'Interpreter', 'latex');
xlim([x_min x_max])
ylim([-1 1])

% Formatting
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 13)
set(gcf, 'Position', [500 100 1000 400])

% Labels
text(h(1), -3.8, 1, '(a)', 'Interpreter', 'latex', 'FontSize', 15)
text(h(2), -3.8, 1, '(b)', 'Interpreter', 'latex', 'FontSize', 15)

% Save options (uncomment if needed)
% saveas(gcf, 'approach_to_terminal.fig')
% saveas(gcf, 'approach_to_terminal.eps', 'epsc')
% saveas(gcf, 'approach_to_terminal.jpg')

%% double wave-period averaging for Figure 9

% from here to the end, analysis using wave-averaging is provided, which
% will be used to draw figure 9.

clc
close all
clear all
load total_wave_filt
N1_W1_data = nylon_332;

h_dep=315;
v_settle=-74.5;

data_length = length(N1_W1_data);
wave_average_2nd = cell(1,data_length);
for kk=1:data_length
    clear test_vv test_zz mean_vv_1st mean_zz_1st mean_vv_2nd mean_zz_2nd
    test_vv = N1_W1_data{1,kk}(:,5);
    test_zz = N1_W1_data{1,kk}(:,3);
    
    num_wave = 13; % waveperiod/sampling frequency.
    for ii=1:length(test_vv)-num_wave
        mean_vv_1st(ii) = mean(test_vv(ii:ii+num_wave-1));
        mean_zz_1st(ii) = mean(test_zz(ii:ii+num_wave-1));
    end
    for ii=1:length(mean_vv_1st)-num_wave
        mean_vv_2nd(ii) = mean(mean_vv_1st(ii:ii+num_wave-1));
        mean_zz_2nd(ii) = mean(mean_zz_1st(ii:ii+num_wave-1));
    end
    wave_average_2nd{1,kk}(:,1) = mean_zz_2nd;
    wave_average_2nd{1,kk}(:,2) = mean_vv_2nd;
end

% drift from analytical model using W1 wave case, Eq. (4.4)

eps = 0.0714;
zp = -(0:0.01:0.8);
kh = 2.95;
k=9.37;
vz_est = eps(1).^2.*cosh(2*(kh.*zp+kh))./(2*(cosh(kh)^2))/(1+(N1_vs_norm)^2);


%% separate double-averaged data for every 0.1z/h

clc
close all
clear v_dep z_dep
velo_save=1;
data_length = length(N1_W1_data);

% take data where (double averaged vertical position) > -0.2 h_p
for kk=1:data_length
    clear velo_save zz_save
    ss = find(wave_average_2nd{1,kk}(:,1)>-0.2*h_dep);
    velo_save = wave_average_2nd{1,kk}(ss,2);
    zz_save = wave_average_2nd{1,kk}(ss,1);
    v_dep{1,kk}(:,1) = velo_save;
    z_dep{1,kk}(:,1) = zz_save;
end

% gather data for bin averaging
num_ini=1;
for kk=1:data_length
    N1_W1_vz_12(num_ini:num_ini+length(z_dep{1,kk})-1,1) = v_dep{1,kk}(:,1);
    N1_W1_zp_12(num_ini:num_ini+length(v_dep{1,kk})-1,1) = z_dep{1,kk}(:,1);
    num_ini= num_ini+length(v_dep{1,kk});
end

% take data where -0.3 h_p < (double averaged vertical position) < -0.2 h_p 

clear v_dep z_dep
for kk=1:data_length
    clear velo_save zz_save
    ss = find(wave_average_2nd{1,kk}(:,1)<-0.2*h_dep & wave_average_2nd{1,kk}(:,1)>-0.3*h_dep);
    velo_save = wave_average_2nd{1,kk}(ss,2);
    zz_save = wave_average_2nd{1,kk}(ss,1);
    v_dep{1,kk}(:,1) = velo_save;
    z_dep{1,kk}(:,1) = zz_save ;
end
% gather data for bin averaging

num_ini = 1;
for kk=1:data_length
    N1_W1_vz_23(num_ini:num_ini+length(v_dep{1,kk})-1,1) = v_dep{1,kk}(:,1);
    N1_W1_zp_23(num_ini:num_ini+length(v_dep{1,kk})-1,1) = z_dep{1,kk}(:,1);
    num_ini= num_ini+length(v_dep{1,kk});
end

% take data where -0.4 h_p < (double averaged vertical position) < -0.3 h_p 

clear v_dep z_dep
for kk=1:data_length
    clear velo_save zz_save
    ss = find(wave_average_2nd{1,kk}(:,1)<-0.3*h_dep & wave_average_2nd{1,kk}(:,1)>-0.4*h_dep);
    velo_save = wave_average_2nd{1,kk}(ss,2);
    zz_save = wave_average_2nd{1,kk}(ss,1);
    v_dep{1,kk}(:,1) = velo_save;
    z_dep{1,kk}(:,1) = zz_save ;
end
% gather data for bin averaging

num_ini = 1;
for kk=1:data_length
    N1_W1_vz_34(num_ini:num_ini+length(v_dep{1,kk})-1,1) = v_dep{1,kk}(:,1);
    N1_W1_zp_34(num_ini:num_ini+length(v_dep{1,kk})-1,1) = z_dep{1,kk}(:,1);
    num_ini= num_ini+length(v_dep{1,kk});
end

% take data where -0.5 h_p < (double averaged vertical position) < -0.4 h_p 

clear v_dep z_dep
for kk=1:data_length
    clear velo_save zz_save
    ss = find(wave_average_2nd{1,kk}<-0.4*h_dep & wave_average_2nd{1,kk}>-0.5*h_dep);
    velo_save = wave_average_2nd{1,kk}(ss,2);
    zz_save = wave_average_2nd{1,kk}(ss,1);
    v_dep{1,kk}(:,1) = velo_save;
    z_dep{1,kk}(:,1) = zz_save ;
end
% gather data for bin averaging

num_ini = 1;
for kk=1:data_length
    N1_W1_vz_45(num_ini:num_ini+length(v_dep{1,kk})-1,1) = v_dep{1,kk}(:,1);
    N1_W1_zp_45(num_ini:num_ini+length(v_dep{1,kk})-1,1) = z_dep{1,kk}(:,1);
    num_ini= num_ini+length(v_dep{1,kk});
end

% take data where -0.6 h_p < (double averaged vertical position) < -0.5 h_p

clear v_dep z_dep
for kk=1:data_length
    clear velo_save zz_save
    ss = find(wave_average_2nd{1,kk}<-0.5*h_dep & wave_average_2nd{1,kk}>-0.6*h_dep);
    velo_save = wave_average_2nd{1,kk}(ss,2);
    zz_save = wave_average_2nd{1,kk}(ss,1);
    v_dep{1,kk}(:,1) = velo_save;
    z_dep{1,kk}(:,1) = zz_save;
end
% gather data for bin averaging

num_ini = 1;
for kk=1:data_length
    N1_W1_vz_56(num_ini:num_ini+length(v_dep{1,kk})-1,1) = v_dep{1,kk}(:,1);
    N1_W1_zp_56(num_ini:num_ini+length(v_dep{1,kk})-1,1) = z_dep{1,kk}(:,1);
    num_ini= num_ini+length(v_dep{1,kk});
end

% take data where -0.7 h_p < (double averaged vertical position) < -0.6 h_p

clear v_dep z_dep
for kk=1:data_length
    clear velo_save zz_save
    ss = find(wave_average_2nd{1,kk}<-0.6*h_dep & wave_average_2nd{1,kk}>-0.7*h_dep);
    velo_save = wave_average_2nd{1,kk}(ss,2);
    zz_save = wave_average_2nd{1,kk}(ss,1);
    v_dep{1,kk}(:,1) = velo_save;
    z_dep{1,kk}(:,1) = zz_save;
   
end

% gather data for bin averaging

num_ini = 1;
for kk=1:data_length
    N1_W1_vz_67(num_ini:num_ini+length(v_dep{1,kk})-1,1) = v_dep{1,kk}(:,1);
    N1_W1_zp_67(num_ini:num_ini+length(v_dep{1,kk})-1,1) = z_dep{1,kk}(:,1);
    num_ini= num_ini+length(v_dep{1,kk});
end

%% bin averaging via bootstrapping
clc
close all
nboots = 2000

mean_N1_W1_zp_12 = mean(bootstrp(nboots,@mean,N1_W1_zp_12));
mean_N1_W1_zp_23 = mean(bootstrp(nboots,@mean,N1_W1_zp_23));
mean_N1_W1_zp_34 = mean(bootstrp(nboots,@mean,N1_W1_zp_34));
mean_N1_W1_zp_45 = mean(bootstrp(nboots,@mean,N1_W1_zp_45));
mean_N1_W1_zp_56 = mean(bootstrp(nboots,@mean,N1_W1_zp_56));
mean_N1_W1_zp_67 = mean(bootstrp(nboots,@mean,N1_W1_zp_67));

ci_N1_W1_zp_12 = bootci(nboots,@mean,N1_W1_zp_12);
ci_N1_W1_zp_23 = bootci(nboots,@mean,N1_W1_zp_23);
ci_N1_W1_zp_34 = bootci(nboots,@mean,N1_W1_zp_34);
ci_N1_W1_zp_45 = bootci(nboots,@mean,N1_W1_zp_45);
ci_N1_W1_zp_56 = bootci(nboots,@mean,N1_W1_zp_56);
ci_N1_W1_zp_67 = bootci(nboots,@mean,N1_W1_zp_67);

mean_N1_W1_zp = [mean_N1_W1_zp_12 mean_N1_W1_zp_23 mean_N1_W1_zp_34 mean_N1_W1_zp_45 mean_N1_W1_zp_56 mean_N1_W1_zp_67];
neg_N1_W1_zp = [ci_N1_W1_zp_12(1) ci_N1_W1_zp_23(1) ci_N1_W1_zp_34(1) ci_N1_W1_zp_45(1) ci_N1_W1_zp_56(1) ci_N1_W1_zp_67(1)];
pos_N1_W1_zp = [ci_N1_W1_zp_12(2) ci_N1_W1_zp_23(2) ci_N1_W1_zp_34(2) ci_N1_W1_zp_45(2) ci_N1_W1_zp_56(2) ci_N1_W1_zp_67(2)];

mean_N1_W1_vz_12 = mean(bootstrp(nboots,@mean,N1_W1_vz_12));
mean_N1_W1_vz_23 = mean(bootstrp(nboots,@mean,N1_W1_vz_23));
mean_N1_W1_vz_34 = mean(bootstrp(nboots,@mean,N1_W1_vz_34));
mean_N1_W1_vz_45 = mean(bootstrp(nboots,@mean,N1_W1_vz_45));
mean_N1_W1_vz_56 = mean(bootstrp(nboots,@mean,N1_W1_vz_56));
mean_N1_W1_vz_67 = mean(bootstrp(nboots,@mean,N1_W1_vz_67));

ci_N1_W1_vz_12 = bootci(nboots,@mean,N1_W1_vz_12);
ci_N1_W1_vz_23 = bootci(nboots,@mean,N1_W1_vz_23);
ci_N1_W1_vz_34 = bootci(nboots,@mean,N1_W1_vz_34);
ci_N1_W1_vz_45 = bootci(nboots,@mean,N1_W1_vz_45);
ci_N1_W1_vz_56 = bootci(nboots,@mean,N1_W1_vz_56);
ci_N1_W1_vz_67 = bootci(nboots,@mean,N1_W1_vz_67);

mean_N1_W1_vz = [mean_N1_W1_vz_12 mean_N1_W1_vz_23 mean_N1_W1_vz_34 mean_N1_W1_vz_45 mean_N1_W1_vz_56 mean_N1_W1_vz_67];
neg_N1_W1_vz = [ci_N1_W1_vz_12(1) ci_N1_W1_vz_23(1) ci_N1_W1_vz_34(1) ci_N1_W1_vz_45(1) ci_N1_W1_vz_56(1) ci_N1_W1_vz_67(1)];
pos_N1_W1_vz = [ci_N1_W1_vz_12(2) ci_N1_W1_vz_23(2) ci_N1_W1_vz_34(2) ci_N1_W1_vz_45(2) ci_N1_W1_vz_56(2) ci_N1_W1_vz_67(2)];


%% dispersion interval = v_s*T_wave/h_dep

clc
close all
clear disp_zz_N1 disp_N1_x pos_x_N1 mean_dep_z_N1 ci_dep_z_N1 


% over measurement area, we consider z'/h_dep < -0.2 (dz_ini_N1) and then
% separate date for every settling distance per single wave period.
dz_ini = 0.0;
dz_ini_N1 = 0.2;
dz_end = 0.8;
L_wave = 670.6;

dz_N1 = 0.157; % v_s*T/h_dep (fraction of settling distance per single wave period)
num_jj = floor((dz_end-dz_ini)/dz_N1);
for kk=1:10
    for jj=1:num_jj
        ii = nylon_332{kk}(:,3)<-(dz_ini+dz_N1*(jj-1))*h_dep & nylon_332{kk}(:,3)>-(dz_ini+dz_N1*jj)*h_dep;
        WA_leng_N1(jj,kk) = sum(ii);
        pos_x_N1(jj,kk) = mean(nylon_332{kk}(ii,2)-nylon_332{kk}(1,2));
        pos_z_N1(jj,kk) = mean(nylon_332{kk}(ii,3));
    end
end

disp_N1_x = std(pos_x_N1,0,2);
disp_N1_z = std(pos_z_N1,0,2);

nboot = 1000;
% mean_z
for ii=1:size(pos_z_N1,1)
    data = pos_z_N1(ii,:);
    m = bootstrp(nboot,@mean,data);
    mean_dep_z_N1(ii,1) = mean(m);
    ci_dep_z_N1(ii,:) = bootci(nboot,@mean,data)-mean(m);
end
 
%% Figure 9: enhanced settling + horizontal dispersion
clc
close all

h(1) = subplot(121);


errorbar((mean_N1_W1_vz)/v_settle,mean_N1_W1_zp/h_dep, (mean_N1_W1_zp-neg_N1_W1_zp)/h_dep,(mean_N1_W1_zp-pos_N1_W1_zp)/h_dep,...
    abs((mean_N1_W1_vz-neg_N1_W1_vz)/v_settle),((mean_N1_W1_vz-pos_N1_W1_vz)/v_settle),'or')
ylim([-0.7 -.1])

hold on
plot([1 1],[-0.8 -0],'-k')
yticks([-0.8:0.1:0]) 
ylim([-0.8 -0])


std_vs = 1;
err_vsm =sqrt(std_vs^2);

plot(1+[1 1]*(err_vsm)/v_settle,[-0.8 -0],'-.k')
plot(1-[1 1]*(err_vsm)/v_settle,[-0.8 -0],'-.k')

plot(1+vz_est,zp,'--r')
xlabel('$ -\overline{\overline{v_z}}/v_s $','interpreter','latex')
ylabel('$ \overline{\overline{z_p}}/h $','interpreter','latex')

h(2) = subplot(122);

errorbar(disp_N1_x/L_wave,mean_dep_z_N1/h_dep, ci_dep_z_N1(:,1)/h_dep, ci_dep_z_N1(:,2)/h_dep,'-or')
hold on

xlabel('$\sigma_{\overline{x_p}}/L$','interpreter','latex')
ylabel('$\overline{z_p}/h$','interpreter','latex')
xlim([0.006 0.01])
ylim([-0.8 -0])


set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf,'position',[100 100 800 400])
text(h(1),0.92,0,'(a)','fontsize',14,'interpreter','latex')
text(h(2),0.005,0,'(b)','fontsize',14,'interpreter','latex')

% optional
% saveas(figure(1),'enhanced_settling_dispersion.fig')
% saveas(figure(1),'enhanced_settling_dispersion.jpg')
% saveas(figure(1),'enhanced_settling_dispersion.eps','epsc')

