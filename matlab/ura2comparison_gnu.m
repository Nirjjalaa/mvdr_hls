clc;
clear;
close all;

% ================= PARAMETERS =================
Nx = 8;
Ny = 8;
N  = Nx * Ny;

d = 0.5;
K = 5000;

% ================= ANGLES =================
az_s  = 40;  el_s  = 0;
az_i1 = 30;  el_i1 = 10;
az_i2 = -60; el_i2 = -25;

% ================= REPRODUCIBILITY =================
rand("seed", 42);

% ================= SIGNALS =================
t = (0:K-1)/K;

s  = cos(2*pi*50*t);
i1 = randn(1,K) + 1j*randn(1,K);
i2 = randn(1,K) + 1j*randn(1,K);

noise = sqrt(0.01) * (randn(N,K) + 1j*randn(N,K));

% ================= STEERING VECTORS =================
a_s  = steering_vec(az_s,  el_s,  Nx, Ny, d);
a_i1 = steering_vec(az_i1, el_i1, Nx, Ny, d);
a_i2 = steering_vec(az_i2, el_i2, Nx, Ny, d);

% ================= DATA MODEL =================
Xtrain = 100*(a_i1*i1) + 100*(a_i2*i2) + noise;
Xdata  = 10*(a_s*s) + 100*(a_i1*i1) + 100*(a_i2*i2) + noise;

% ================= COVARIANCE MATRIX =================
R = (Xtrain * Xtrain') / K;

delta = 1e-4 * trace(R) / N;
R_reg = R + delta * eye(N);

% ================= TIMING ARRAYS =================
NUM_TRIALS = 50;

t_builtin = zeros(1, NUM_TRIALS);
t_manual  = zeros(1, NUM_TRIALS);

% =========================================================
% ================= BUILT-IN MVDR ==========================
% =========================================================

for tr = 1:NUM_TRIALS
    tic;

    w_builtin = R_reg \ a_s;   % fast solver (Octave-supported)
    w_builtin = w_builtin / (a_s' * w_builtin);

    t_builtin(tr) = toc;
end

% =========================================================
% ================= MANUAL MVDR ============================
% =========================================================

for tr = 1:NUM_TRIALS
    tic;

    [Q, Rq] = givens_qr(R_reg);

    y = Q' * a_s;

    w_manual = back_substitution(Rq, y);

    w_manual = w_manual / (a_s' * w_manual);

    t_manual(tr) = toc;
end

% =========================================================
% ================= TIMING RESULTS =========================
% =========================================================

fprintf("\n===== TIMING RESULTS (%d trials) =====\n", NUM_TRIALS);

fprintf("Built-in  mean: %.4f ms | std: %.4f ms\n", ...
    mean(t_builtin)*1e3, std(t_builtin)*1e3);

fprintf("Manual    mean: %.4f ms | std: %.4f ms\n", ...
    mean(t_manual)*1e3, std(t_manual)*1e3);

fprintf("Overhead factor (Manual/Built-in): %.2fx\n", ...
    mean(t_manual)/mean(t_builtin));

% =========================================================
% ================= GAIN CHECK ============================
% =========================================================

g_s    = abs(w_manual' * a_s);
g_i1   = abs(w_manual' * a_i1);
g_i2   = abs(w_manual' * a_i2);

g_s_b  = abs(w_builtin' * a_s);
g_i1_b = abs(w_builtin' * a_i1);
g_i2_b = abs(w_builtin' * a_i2);

fprintf("\n===== MVDR GAINS (Manual vs Built-in) =====\n");

fprintf("Target Gain: %.4f | %.4f\n", g_s, g_s_b);
fprintf("I1 Gain:     %.6f | %.6f\n", g_i1, g_i1_b);
fprintf("I2 Gain:     %.6f | %.6f\n", g_i2, g_i2_b);

fprintf("\nSuppression (dB):\n");

fprintf("I1 rejection: %.2f dB | %.2f dB\n", ...
    20*log10(g_i1/g_s), 20*log10(g_i1_b/g_s_b));

fprintf("I2 rejection: %.2f dB | %.2f dB\n", ...
    20*log10(g_i2/g_s), 20*log10(g_i2_b/g_s_b));

% =========================================================
% ================= 2D BEAM PATTERN ========================
% =========================================================

az_scan = -90:1:90;
el_scan = -90:1:90;

P = zeros(length(el_scan), length(az_scan));

for ii = 1:length(el_scan)
    for jj = 1:length(az_scan)

        a = steering_vec(az_scan(jj), el_scan(ii), Nx, Ny, d);

        P(ii,jj) = abs(w_manual' * a);

    end
end

% Normalize
P = P / max(P(:));
PdB = 20*log10(P + 1e-12);

% =========================================================
% ================= PLOT ================================
% =========================================================

figure;
imagesc(az_scan, el_scan, PdB);
axis xy;
colorbar;
caxis([-40 0]);   % Octave-safe replacement for clim()

xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
title('URA MVDR Beam Pattern (Manual)');

hold on;

% ================= MARKERS =================
plot(az_s,  el_s,  'g^', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(az_i1, el_i1, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(az_i2, el_i2, 'mv', 'MarkerSize', 10, 'MarkerFaceColor', 'm');

legend('Target','Interference 1','Interference 2');

% =========================================================
% ================= AZIMUTH CUT ============================
% =========================================================

theta_scan = -90:0.5:90;

beam_az_m = zeros(size(theta_scan));
beam_az_b = zeros(size(theta_scan));

for k = 1:length(theta_scan)

    a = steering_vec(theta_scan(k), 0, Nx, Ny, d);

    beam_az_m(k) = abs(w_manual'  * a);
    beam_az_b(k) = abs(w_builtin' * a);

end

% Normalize
beam_az_m = beam_az_m / max(beam_az_m);
beam_az_b = beam_az_b / max(beam_az_b);

beam_az_mdB = 20*log10(beam_az_m + 1e-12);
beam_az_bdB = 20*log10(beam_az_b + 1e-12);

% =========================================================
% ================= PLOT ================================
% =========================================================

figure;

plot(theta_scan, beam_az_mdB, 'b-', 'LineWidth', 2);
hold on;
plot(theta_scan, beam_az_bdB, 'r--', 'LineWidth', 2);

grid on;

xlabel('Azimuth (deg)');
ylabel('Gain (dB)');
title('URA MVDR - Azimuth Cut (Elevation = 0°)');

ylim([-80 5]);

legend('Manual (Givens QR)', 'Built-in (\)');

% =========================================================
% ======= Octave-safe replacements for xline() ===========
% =========================================================

yl = ylim;

plot([az_s az_s], yl, 'g--', 'LineWidth', 1.5);
plot([az_i1 az_i1], yl, 'm--', 'LineWidth', 1.5);
plot([az_i2 az_i2], yl, 'r--', 'LineWidth', 1.5);

text(az_s+2,  -10, 'Target');
text(az_i1+2, -20, 'I1');
text(az_i2+2, -30, 'I2');

function a = steering_vec(az, el, Nx, Ny, d)

    az = az * pi/180;
    el = el * pi/180;

    u = sin(az) * cos(el);
    v = sin(el);

    a = zeros(Nx*Ny, 1);

    k = 1;

    for y = 0:Ny-1
        for x = 0:Nx-1

            phase = 2*pi*d*(x*u + y*v);

            a(k) = exp(1j * phase);
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

            if abs(R(i,j)) > 1e-12

                a = R(i-1,j);
                b = R(i,j);

                r = sqrt(abs(a)^2 + abs(b)^2);

                c = a / r;
                s = -b / r;

                G = [c -s;
                     s  c];

                R([i-1 i], j:n) = G * R([i-1 i], j:n);
                Q(:, [i-1 i])   = Q(:, [i-1 i]) * G';

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
            x(i) = x(i) - R(i,j) * x(j);
        end

        x(i) = x(i) / R(i,i);

    end

end


