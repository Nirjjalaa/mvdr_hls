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
az_s  = 40;
el_s  = 0;

az_i1 = 30;
el_i1 = 10;

az_i2 = -60;
el_i2 = -25;

% ================= SIGNALS =================
t = (0:K-1)/K;

s  = cos(2*pi*50*t);

i1 = randn(1,K) + 1j*randn(1,K);
i2 = randn(1,K) + 1j*randn(1,K);

noise = sqrt(0.5)*(randn(N,K) + 1j*randn(N,K));

% ================= STEERING VECTORS =================
a_s  = steering_vec(az_s,el_s,Nx,Ny,d);
a_i1 = steering_vec(az_i1,el_i1,Nx,Ny,d);
a_i2 = steering_vec(az_i2,el_i2,Nx,Ny,d);

% ================= DATA =================
Xtrain = 5*(a_i1*i1) + 5*(a_i2*i2) + noise;

Xdata = 10*(a_s*s) + ...
        5*(a_i1*i1) + ...
        5*(a_i2*i2) + ...
        noise;

% ================= COVARIANCE =================
R = (Xtrain*Xtrain')/K;

delta = 0.01*trace(R)/N;
R = R + delta*eye(N);

% ================= MANUAL QR =================
[Q,Ru] = givens_qr(R);

y = Q' * a_s;

w_manual = back_substitution(Ru,y);

w_manual = w_manual/(a_s'*w_manual);

% ================= GAIN CHECK =================

g_target = abs(w_manual' * a_s);
g_i1 = abs(w_manual' * a_i1);
g_i2 = abs(w_manual' * a_i2);

disp("===== MVDR GAINS =====");

fprintf("Target Gain (should be ~1): %.4f\n", g_target);
fprintf("Interference 1 Gain: %.6f\n", g_i1);
fprintf("Interference 2 Gain: %.6f\n", g_i2);

fprintf("\nSuppression (dB):\n");
fprintf("I1 rejection: %.2f dB\n",20*log10(g_i1/g_target));
fprintf("I2 rejection: %.2f dB\n",20*log10(g_i2/g_target));

% ================= 2D BEAM PATTERN =================

az = -90:2:90;
el = -90:2:90;

P = zeros(length(el),length(az));

for i = 1:length(el)

    for j = 1:length(az)

        a = steering_vec(az(j),el(i),Nx,Ny,d);

        P(i,j) = abs(w_manual' * a);

    end

end

P = P/max(P(:));

PdB = 20*log10(P + 1e-12);

figure;
imagesc(az,el,PdB);

axis xy;
colorbar;
caxis([-40 0]);

xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
title('Manual URA MVDR Beam Pattern');

drawnow;

% ================= AZIMUTH CUT =================

theta_scan = -80:0.5:80;
beam_az = zeros(size(theta_scan));

for k = 1:length(theta_scan)

    a = steering_vec(theta_scan(k), 0, Nx, Ny, d);
    beam_az(k) = abs(w_manual' * a);

end

beam_az = beam_az / max(beam_az);
beam_az_dB = 20*log10(beam_az + 1e-12);

figure;
plot(theta_scan, beam_az_dB, 'b', 'LineWidth', 2);
grid on;
hold on;

xlabel('Azimuth (deg)');
ylabel('Gain (dB)');
title('URA MVDR - Azimuth Cut (Elevation = 0°)');
ylim([-60 0]);

% ---------- Target & Interference Markers ----------

yl = ylim;

plot([az_s az_s], yl, 'g--', 'LineWidth', 1.5);
plot([az_i1 az_i1], yl, 'r--', 'LineWidth', 1.5);
plot([az_i2 az_i2], yl, 'm--', 'LineWidth', 1.5);

text(az_s+2, -5,  'Target');
text(az_i1+2,-10, 'Int-1');
text(az_i2+2,-15, 'Int-2');

legend('MVDR Pattern','Location','SouthWest');

drawnow;

% ================= ELEVATION CUT =================

el_scan = -80:0.5:80;
beam_el = zeros(size(el_scan));

for k = 1:length(el_scan)

    a = steering_vec(0, el_scan(k), Nx, Ny, d);
    beam_el(k) = abs(w_manual' * a);

end

beam_el = beam_el / max(beam_el);
beam_el_dB = 20*log10(beam_el + 1e-12);

figure;
plot(el_scan, beam_el_dB, 'b', 'LineWidth', 2);
grid on;
hold on;

xlabel('Elevation (deg)');
ylabel('Gain (dB)');
title('URA MVDR - Elevation Cut (Azimuth = 0°)');
ylim([-60 0]);

% ---------- Target & Interference Markers ----------

yl = ylim;

plot([el_s el_s], yl, 'g--', 'LineWidth', 1.5);
plot([el_i1 el_i1], yl, 'r--', 'LineWidth', 1.5);
plot([el_i2 el_i2], yl, 'm--', 'LineWidth', 1.5);

text(el_s+2, -5,  'Target');
text(el_i1+2,-10, 'Int-1');
text(el_i2+2,-15, 'Int-2');

legend('MVDR Pattern','Location','SouthWest');

drawnow;

% =====================================================
% Steering Vector Function
% =====================================================

function a = steering_vec(az, el, Nx, Ny, d)

az = az*pi/180;
el = el*pi/180;

u = sin(az)*cos(el);
v = sin(el);

a = zeros(Nx*Ny,1);

k = 1;

for y = 0:Ny-1
    for x = 0:Nx-1

        phase = 2*pi*d*(x*u + y*v);

        a(k) = exp(1j*phase);

        k = k + 1;

    end
end

a = a/sqrt(Nx*Ny);

end


% =====================================================
% Manual QR using Givens Rotations
% =====================================================

function [Q,R] = givens_qr(A)

[m,n] = size(A);

Q = eye(m);
R = A;

for j = 1:n

    for i = m:-1:(j+1)

        if abs(R(i,j)) > 1e-12

            a = R(i-1,j);
            b = R(i,j);

            r = sqrt(abs(a)^2 + abs(b)^2);

            c = a/r;
            s = -b/r;

            G = [c -s;
                 s  c];

            R([i-1 i],j:n) = G * R([i-1 i],j:n);
            Q(:,[i-1 i]) = Q(:,[i-1 i]) * G';

        end

    end

end

end


% =====================================================
% Back Substitution
% =====================================================

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
