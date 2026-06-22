clc;
clear;
close all;

% ================= PARAMETERS =================
Nx = 8; Ny = 8; N = Nx * Ny;
d = 0.5;
K = 5000;

% ================= ANGLES =================
az_s  = 40;  el_s  = 0;
az_i1 = 30;  el_i1 = 10;
az_i2 = -60; el_i2 = -25;

% ================= SIGNALS =================
rng(42);  % Reproducibility
t = (0:K-1)/K;
s  = cos(2*pi*50*t);
i1 = randn(1,K) + 1j*randn(1,K);
i2 = randn(1,K) + 1j*randn(1,K);
noise = sqrt(0.01)*(randn(N,K) + 1j*randn(N,K));  % Low noise floor

a_s  = steering_vec(az_s,  el_s,  Nx, Ny, d);
a_i1 = steering_vec(az_i1, el_i1, Nx, Ny, d);
a_i2 = steering_vec(az_i2, el_i2, Nx, Ny, d);

% High INR for guaranteed deep nulls
Xtrain = 100*(a_i1*i1) + 100*(a_i2*i2) + noise;
Xdata  = 10*(a_s*s) + 100*(a_i1*i1) + 100*(a_i2*i2) + noise;

% ================= COVARIANCE =================
R = (Xtrain * Xtrain') / K;
delta = 1e-4 * trace(R) / N;   % Small diagonal loading
R_reg = R + delta * eye(N);

% ================= TIMING: BUILT-IN vs MANUAL =================
NUM_TRIALS = 50;

% --- Built-in (MATLAB backslash / mldivide) ---
t_builtin = zeros(1, NUM_TRIALS);
for tr = 1:NUM_TRIALS
    tic;
    w_builtin = R_reg \ a_s;
    w_builtin = w_builtin / (a_s' * w_builtin);
    t_builtin(tr) = toc;
end

% --- Manual (Givens QR + Back Substitution) ---
t_manual = zeros(1, NUM_TRIALS);
for tr = 1:NUM_TRIALS
    tic;
    [Q, Ru] = givens_qr(R_reg);
    y = Q' * a_s;
    w_manual = back_substitution(Ru, y);
    w_manual = w_manual / (a_s' * w_manual);
    t_manual(tr) = toc;
end

fprintf("\n===== TIMING RESULTS (%d trials) =====\n", NUM_TRIALS);
fprintf("Built-in  mean: %.4f ms | std: %.4f ms\n", mean(t_builtin)*1e3, std(t_builtin)*1e3);
fprintf("Manual    mean: %.4f ms | std: %.4f ms\n", mean(t_manual)*1e3,  std(t_manual)*1e3);
fprintf("Overhead factor: %.2fx\n", mean(t_manual)/mean(t_builtin));

% ================= GAIN CHECK =================
g_s    = abs(w_manual' * a_s);
g_i1   = abs(w_manual' * a_i1);
g_i2   = abs(w_manual' * a_i2);
g_s_b  = abs(w_builtin' * a_s);
g_i1_b = abs(w_builtin' * a_i1);
g_i2_b = abs(w_builtin' * a_i2);

fprintf("\n===== MVDR GAINS =====\n");
fprintf("                    Manual     Built-in\n");
fprintf("Target Gain (~1): %.4f     %.4f\n",  g_s,  g_s_b);
fprintf("I1 Gain:          %.6f   %.6f\n",  g_i1, g_i1_b);
fprintf("I2 Gain:          %.6f   %.6f\n",  g_i2, g_i2_b);
fprintf("\nSuppression (dB):\n");
fprintf("I1 rejection: %.2f dB (manual) | %.2f dB (built-in)\n", ...
    20*log10(g_i1/g_s), 20*log10(g_i1_b/g_s_b));
fprintf("I2 rejection: %.2f dB (manual) | %.2f dB (built-in)\n", ...
    20*log10(g_i2/g_s), 20*log10(g_i2_b/g_s_b));

% ================= FIGURE 1: TIMING COMPARISON BAR CHART =================
figure('Name','Timing Comparison','Position',[100 100 700 450]);
methods = {'Built-in (\backslash)', 'Manual (Givens QR)'};
means_ms = [mean(t_builtin), mean(t_manual)] * 1e3;
stds_ms  = [std(t_builtin),  std(t_manual)]  * 1e3;

b = bar(means_ms, 0.5, 'FaceColor','flat');
b.CData(1,:) = [0.2 0.6 0.9];
b.CData(2,:) = [0.9 0.4 0.2];
hold on;
errorbar(1:2, means_ms, stds_ms, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
set(gca, 'XTickLabel', methods, 'FontSize', 12);
ylabel('Time (ms)', 'FontSize', 12);
title(sprintf('Execution Time: %d Trials', NUM_TRIALS), 'FontSize', 13);
grid on;
for k = 1:2
    text(k, means_ms(k) + stds_ms(k) + 0.01*max(means_ms), ...
        sprintf('%.3f ms', means_ms(k)), ...
        'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
end

% ================= FIGURE 2: TIMING TRIAL-BY-TRIAL COMPARISON =================
figure('Name','Trial-by-Trial Timing','Position',[100 100 800 400]);
plot(1:NUM_TRIALS, t_builtin*1e3, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2, ...
    'DisplayName', 'Built-in');
hold on;
plot(1:NUM_TRIALS, t_manual*1e3, 'r-s', 'MarkerSize', 3, 'LineWidth', 1.2, ...
    'DisplayName', 'Manual (Givens QR)');
xlabel('Trial Number', 'FontSize', 12);
ylabel('Time (ms)', 'FontSize', 12);
title('Execution Time per Trial', 'FontSize', 13);
legend('FontSize', 11);
grid on;

% ================= FIGURE 3: 2D BEAM PATTERN HEATMAP (MANUAL) =================
az_scan = -90:1:90;
el_scan = -90:1:90;
P = zeros(length(el_scan), length(az_scan));

for ii = 1:length(el_scan)
    for jj = 1:length(az_scan)
        a = steering_vec(az_scan(jj), el_scan(ii), Nx, Ny, d);
        P(ii,jj) = abs(w_manual' * a);
    end
end
P  = P / max(P(:));
PdB = 20*log10(P + 1e-12);

figure('Name','2D Beam Pattern','Position',[100 100 800 550]);
imagesc(az_scan, el_scan, PdB);
axis xy; colorbar; clim([-40 0]);
hold on;
plot(az_s,  el_s,  'g^', 'MarkerSize', 12, 'MarkerFaceColor','g',  'DisplayName','Target');
plot(az_i1, el_i1, 'rv', 'MarkerSize', 12, 'MarkerFaceColor','r',  'DisplayName','I1');
plot(az_i2, el_i2, 'rv', 'MarkerSize', 12, 'MarkerFaceColor','m',  'DisplayName','I2');
legend('Location','northeast','FontSize',10);
xlabel('Azimuth (deg)','FontSize',12);
ylabel('Elevation (deg)','FontSize',12);
title('URA MVDR 2D Beam Pattern (Manual)','FontSize',13);

% ================= FIGURE 4: AZIMUTH CUT (el = 0) =================
theta_scan = -90:0.5:90;
beam_az_m = zeros(size(theta_scan));
beam_az_b = zeros(size(theta_scan));

for k = 1:length(theta_scan)
    a = steering_vec(theta_scan(k), 0, Nx, Ny, d);
    beam_az_m(k) = abs(w_manual'  * a);
    beam_az_b(k) = abs(w_builtin' * a);
end

beam_az_m  = beam_az_m / max(beam_az_m);
beam_az_b  = beam_az_b / max(beam_az_b);
beam_az_mdB = 20*log10(beam_az_m + 1e-12);
beam_az_bdB = 20*log10(beam_az_b + 1e-12);

figure('Name','Azimuth Cut','Position',[100 100 850 420]);
plot(theta_scan, beam_az_mdB, 'b-',  'LineWidth', 2, 'DisplayName', 'Manual (Givens QR)');
hold on;
plot(theta_scan, beam_az_bdB, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Built-in (\backslash)');
xline(az_s,  'g-',  'LineWidth', 1.5, 'Label', sprintf('Target %d°', az_s));
xline(az_i1, 'm--', 'LineWidth', 1.5, 'Label', sprintf('I1 %d°', az_i1));
xline(az_i2, 'r--', 'LineWidth', 1.5, 'Label', sprintf('I2 %d°', az_i2));
grid on;
xlabel('Azimuth (deg)', 'FontSize', 12);
ylabel('Gain (dB)', 'FontSize', 12);
title('URA MVDR – Azimuth Cut (Elevation = 0°)', 'FontSize', 13);
ylim([-80 5]);
legend('Location','southwest','FontSize',10);

% ================= FIGURE 5: ELEVATION CUT (az = 0) =================
el_cut = -90:0.5:90;
beam_el_m = zeros(size(el_cut));
beam_el_b = zeros(size(el_cut));

for k = 1:length(el_cut)
    a = steering_vec(0, el_cut(k), Nx, Ny, d);
    beam_el_m(k) = abs(w_manual'  * a);
    beam_el_b(k) = abs(w_builtin' * a);
end

beam_el_m   = beam_el_m / max(beam_el_m);
beam_el_b   = beam_el_b / max(beam_el_b);
beam_el_mdB = 20*log10(beam_el_m + 1e-12);
beam_el_bdB = 20*log10(beam_el_b + 1e-12);

figure('Name','Elevation Cut','Position',[100 100 850 420]);
plot(el_cut, beam_el_mdB, 'b-',  'LineWidth', 2, 'DisplayName', 'Manual (Givens QR)');
hold on;
plot(el_cut, beam_el_bdB, 'r--', 'LineWidth', 1.8, 'DisplayName', 'Built-in (\backslash)');
xline(el_s,   'g-',  'LineWidth', 1.5, 'Label', sprintf('Target el=%d°', el_s));
xline(el_i1,  'm--', 'LineWidth', 1.5, 'Label', sprintf('I1 el=%d°', el_i1));
xline(el_i2,  'r--', 'LineWidth', 1.5, 'Label', sprintf('I2 el=%d°', el_i2));
grid on;
xlabel('Elevation (deg)', 'FontSize', 12);
ylabel('Gain (dB)', 'FontSize', 12);
title('URA MVDR – Elevation Cut (Azimuth = 0°)', 'FontSize', 13);
ylim([-80 5]);
legend('Location','southwest','FontSize',10);

% ================= FIGURE 6: OUTPUT SINR COMPARISON =================
sinr_manual  = 10*log10(abs(w_manual'  * a_s)^2 / ...
    (abs(w_manual'  * a_i1)^2 + abs(w_manual'  * a_i2)^2 + 0.01));
sinr_builtin = 10*log10(abs(w_builtin' * a_s)^2 / ...
    (abs(w_builtin' * a_i1)^2 + abs(w_builtin' * a_i2)^2 + 0.01));

figure('Name','SINR Comparison','Position',[100 100 500 400]);
bar_data = [sinr_manual, sinr_builtin];
b2 = bar(bar_data, 0.45, 'FaceColor','flat');
b2.CData(1,:) = [0.2 0.6 0.9];
b2.CData(2,:) = [0.9 0.4 0.2];
set(gca,'XTickLabel',{'Manual (Givens QR)','Built-in (\backslash)'},'FontSize',12);
ylabel('SINR (dB)','FontSize',12);
title('Output SINR Comparison','FontSize',13);
grid on;
for k = 1:2
    text(k, bar_data(k)+0.5, sprintf('%.1f dB', bar_data(k)), ...
        'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
end

% ================= HELPER FUNCTIONS =================
function a = steering_vec(az, el, Nx, Ny, d)
    az_r = deg2rad(az);
    el_r = deg2rad(el);
    u = sin(az_r) * cos(el_r);
    v = sin(el_r);
    a = zeros(Nx*Ny, 1);
    k = 1;
    for iy = 0:Ny-1
        for ix = 0:Nx-1
            a(k) = exp(1j * 2*pi*d * (ix*u + iy*v));
            k = k + 1;
        end
    end
    a = a / sqrt(Nx*Ny);
end

function [Q, R] = givens_qr(A)
    [m, n] = size(A);
    Q = eye(m);
    R = A;
    for j = 1:n
        for i = m:-1:(j+1)
            if abs(R(i,j)) > 1e-15
                a = R(i-1,j); b = R(i,j);
                r = sqrt(a^2 + b^2);
                c = a/r;  s_val = -b/r;
                G = [c -s_val; s_val c];
                R([i-1 i], j:n)  = G  * R([i-1 i], j:n);
                Q(:, [i-1 i])    = Q(:, [i-1 i]) * G';
            end
        end
    end
end

function x = back_substitution(R, b)
    n = length(b);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = b(i);
        for j = i+1:n
            x(i) = x(i) - R(i,j)*x(j);
        end
        x(i) = x(i) / R(i,i);
    end
end
