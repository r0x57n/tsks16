clear; close all; clc;

N = 8192;
n = 0:N-1;

omega1T_coherent = 0.25*pi;
omega1T_noncoherent = (0.25 + 1/N)*pi;
delta_omegaT = 2*pi/N;

p = 50;
amp = 0.01;

window_rect = rectwin(N)';
window_bmh = blackmanharris(N)';

% p represents the number of frequency bins between the two...
% f2 - f1 = p/(T*N);

% Rectangle window
x_coh = sin(omega1T_coherent*n) + amp*sin((omega1T_coherent + p*delta_omegaT)*n);
X_coh = fft(x_coh.*window_rect, N);
X_coh_dB = 20*log10(abs(X_coh)/max(abs(X_coh)));

x_non = sin(omega1T_noncoherent*n) + amp*sin((omega1T_noncoherent + p*delta_omegaT)*n);
X_non = fft(x_non.*window_rect, N);
X_non_dB = 20*log10(abs(X_non)/max(abs(X_non)));

subplot(2,2,1);
plot(X_coh_dB); hold on;
plot(X_non_dB); hold on;
xlim([0 N]);
title("Rectangle window");
ylabel("A (dB)")
xlabel("F (Hz)")

% Blackman Harris window
x_coh = sin(omega1T_coherent*n) + amp*sin((omega1T_coherent + p*delta_omegaT)*n);
X_coh = fft(x_coh.*window_bmh, N);
X_coh_dB = 20*log10(abs(X_coh)/max(abs(X_coh)));

x_non = sin(omega1T_noncoherent*n) + amp*sin((omega1T_noncoherent + p*delta_omegaT)*n);
X_non = fft(x_non.*window_bmh, N);
X_non_dB = 20*log10(abs(X_non)/max(abs(X_non)));

subplot(2,2,2);
plot(X_coh_dB); hold on;
plot(X_non_dB); hold on;
xlim([0 N]);
title("BMH window");
ylabel("A (dB)")
xlabel("F (Hz)")


amp = 0.001;

% Rectangle window
x_coh = sin(omega1T_coherent*n) + amp*sin((omega1T_coherent + p*delta_omegaT)*n);
X_coh = fft(x_coh.*window_rect, N);
X_coh_dB = 20*log10(abs(X_coh)/max(abs(X_coh)));

x_non = sin(omega1T_noncoherent*n) + amp*sin((omega1T_noncoherent + p*delta_omegaT)*n);
X_non = fft(x_non.*window_rect, N);
X_non_dB = 20*log10(abs(X_non)/max(abs(X_non)));

subplot(2,2,3);
plot(X_coh_dB); hold on;
plot(X_non_dB); hold on;
xlim([0 N]);
title("Rectangle window A small");
ylabel("A (dB)")
xlabel("F (Hz)")

% Blackman Harris window
x_coh = sin(omega1T_coherent*n) + amp*sin((omega1T_coherent + p*delta_omegaT)*n);
X_coh = fft(x_coh.*window_bmh, N);
X_coh_dB = 20*log10(abs(X_coh)/max(abs(X_coh)));

x_non = sin(omega1T_noncoherent*n) + amp*sin((omega1T_noncoherent + p*delta_omegaT)*n);
X_non = fft(x_non.*window_bmh, N);
X_non_dB = 20*log10(abs(X_non)/max(abs(X_non)));

subplot(2,2,4);
plot(X_coh_dB); hold on;
plot(X_non_dB); hold on;
xlim([0 N]);
title("BMH window A small");
ylabel("A (dB)")
xlabel("F (Hz)")
