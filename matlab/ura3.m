clc;
clear;
close all;

% ================= PARAMETERS =================
Nx = 8;
Ny = 8;
N = Nx * Ny;

d = 0.5;

K = 5000;

% ================= ANGLES =================
az_s = 40;  el_s = 0;
az_i1 = 30; el_i1 = 10;
az_i2 = -60; el_i2 = -25;

% ================= SIGNALS =================
t = (0:K-1)/K;
s = exp(1j*2*pi*50*t);
i1 = randn(1,K) + 1j*randn(1,K);
i2 = randn(1,K) + 1j*randn(1,K);

noise = sqrt(0.5)*(randn(N,K) + 1j*randn(N,K));

a_s  = steering_vec(az_s,  el_s,  Nx, Ny, d);
a_i1 = steering_vec(az_i1, el_i1, Nx, Ny, d);
a_i2 = steering_vec(az_i2, el_i2, Nx, Ny, d);
Xtrain = 5*(a_i1*i1) + 5*(a_i2*i2) + noise;
Xdata  = 10*(a_s*s)   + 5*(a_i1*i1) + 5*(a_i2*i2) + noise;

R = (Xtrain * Xtrain') / K;

delta = 0.01 * trace(R)/N;
R = R + delta*eye(N);

[Q,Ru] = givens_qr(R);
disp('----------------------------------');
disp('Complex QR Verification');

fprintf('||Q''Q-I|| = %e\n', norm(Q'*Q-eye(size(Q))));
fprintf('||QR-R||  = %e\n', norm(Q*Ru-R));
fprintf('Lower triangular error = %e\n', ...
    norm(tril(Ru,-1)));

disp('----------------------------------');

y = Q' * a_s;

w_manual = back_substitution(Ru, y);
w_manual = w_manual / (a_s' * w_manual);
% ================= GAIN CHECK =================

g_target = abs(w_manual' * a_s);
g_i1 = abs(w_manual' * a_i1);
g_i2 = abs(w_manual' * a_i2);

disp("===== MVDR GAINS =====");

fprintf("Target Gain (should be ~1): %.4f\n", g_target);
fprintf("Interference 1 Gain: %.6f\n", g_i1);
fprintf("Interference 2 Gain: %.6f\n", g_i2);

% Optional suppression in dB
fprintf("\nSuppression (dB):\n");
fprintf("I1 rejection: %.2f dB\n", 20*log10(g_i1/g_target));
fprintf("I2 rejection: %.2f dB\n", 20*log10(g_i2/g_target));
az = -90:0.5:90;
el = -90:0.5:90;

P = zeros(length(el), length(az));

for i = 1:length(el)
    for j = 1:length(az)

        a = steering_vec(az(j), el(i), Nx, Ny, d);
        P(i,j) = abs(w_manual' * a);

    end
end

P = P / max(P(:));
PdB = 20*log10(P + 1e-12);

figure;
imagesc(az, el, PdB);
axis xy;
colorbar;
caxis([-40 0]);

xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
title('Manual URA MVDR (No Toolbox)');


function a = steering_vec(az, el, Nx, Ny, d)

az = deg2rad(az);
el = deg2rad(el);

a = zeros(Nx*Ny,1);

% direction cosines (IMPORTANT FIX)
u = sin(az) * cos(el);
v = sin(el);

k = 1;

for y = 0:Ny-1
    for x = 0:Nx-1

        % CLEAN URA PHYSICS MODEL
        phase = 2*pi*d*(x*u + y*v);

        a(k) = exp(1j*phase);
        k = k + 1;

    end
end

a = a / sqrt(Nx*Ny);

end

function [Q,R] = givens_qr(A)

[m,n] = size(A);

Q = complex(eye(m));
R = complex(A);

for j = 1:n

    for i = m:-1:(j+1)

        if abs(R(i,j)) > 1e-12

            a = R(i-1,j);
            b = R(i,j);

            % Compute norm
            r = sqrt(abs(a)^2 + abs(b)^2);

            if r == 0
                continue;
            end

            % Complex Givens coefficients
            c = abs(a)/r;

            if abs(a) < 1e-12
                alpha = 1;
            else
                alpha = a/abs(a);      % phase of a
            end

            s = alpha * conj(b) / r;

            % Complex Givens rotation
            G = [ c          s;
                 -conj(s)    c ];

            % Apply rotation
            R([i-1 i],j:n) = G * R([i-1 i],j:n);

            % Accumulate Q
            Q(:,[i-1 i]) = Q(:,[i-1 i]) * G';

        end

    end

end

end

function x = back_substitution(R,b)

n = length(b);
x = zeros(n,1);

for i = n:-1:1
    x(i) = b(i);

    for j = i+1:n
        x(i) = x(i) - R(i,j)*x(j);
    end

    x(i) = x(i)/R(i,i);
end

end
theta_scan = -80:0.5:80;
beam_az = zeros(size(theta_scan));

for k = 1:length(theta_scan)

    a = steering_vec(theta_scan(k), 0, Nx, Ny, d);
    beam_az(k) = abs(w_manual' * a);

end

beam_az = beam_az / max(beam_az);
beam_az_dB = 20*log10(beam_az + 1e-12);

figure;
plot(theta_scan, beam_az_dB, 'LineWidth', 2);
grid on;
xlabel('Azimuth (deg)');
ylabel('Gain (dB)');
title('URA MVDR - Azimuth Cut (Elevation = 0°)');
ylim([-60 0]);

hold on;
xline(az_s,'g--','Target');
xline(az_i1,'r--','Interference');
xline(az_i2,'r--','Interference');
el_scan = -80:0.5:80;
beam_el = zeros(size(el_scan));

for k = 1:length(el_scan)

    a = steering_vec(0, el_scan(k), Nx, Ny, d);
    beam_el(k) = abs(w_manual' * a);

end

beam_el = beam_el / max(beam_el);
beam_el_dB = 20*log10(beam_el + 1e-12);

figure;
plot(el_scan, beam_el_dB, 'LineWidth', 2);
grid on;
xlabel('Elevation (deg)');
ylabel('Gain (dB)');
title('URA MVDR - Elevation Cut (Azimuth = 0°)');
ylim([-60 0]);

disp('=== Intermediate Results ===');

fprintf('Steering vector size : %d x %d\n', size(a_s));
fprintf('Covariance matrix R : %d x %d\n', size(R));
fprintf('Q matrix            : %d x %d\n', size(Q));
fprintf('Upper matrix Ru     : %d x %d\n', size(Ru));
fprintf('Weight vector       : %d x %d\n', size(w_manual));

fprintf('\nNorm of weight vector = %.6f\n', norm(w_manual));
fprintf('Condition number of R = %.2f\n', cond(R));

save('golden_reference.mat', ...
    'a_s','R','Ru','Q','y','w_manual');
disp('First 15 elements of steering vector:');
disp(a_s(1:15));

disp('Top-left 5x5 of R:');
disp(R(1:5,1:5));

disp('Top-left 5x5 of Q:');
disp(Q(1:5,1:5));

disp('Top-left 5x5 of Ru:');
disp(Ru(1:5,1:5));

disp('First 15 weights:');
disp(w_manual(1:15));
