clear; close all; clc;

Fs = 1e3;
Ts = 1/Fs;
t = 0:Ts:1;

f1 = 45.*rand();
f2 = 45.*rand();
x = 1.5*cos(2*pi*f1*t) + 2*cos(2*pi*f2*t);

%% lab2_1 test
plot(x); hold on;
title("unity-gain AGC");
x_agc = AGCunity(x);
plot(x_agc); hold on;
line([0, 1e3], [1, 1], "color", "black")
line([0, 1e3], [-1, -1], "color", "black")
legend("Original", "AGC");
xlim([0 1e3]);
grid on;

%% lab2_2 test
figure;
xqs = zeros(5, length(x));
xqs(1,:) = quant(x_agc, 1);
xqs(2,:) = quant(x_agc, 2);
xqs(3,:) = quant(x_agc, 4);
xqs(4,:) = quant(x_agc, 8);
xqs(5,:) = quant(x_agc, 16);
plot(x_agc); hold on;
title("Quantization");

for i = 0:4
    subplot(5,1,i+1);
    plot(x_agc); hold on;
    plot(xqs(i+1,:)); hold on;
    title("B="+2^i);
    xlim([0 1e3]);
end

%% lab2_3
N = 2^20;
n = 0:N-1;
xes = zeros(4,N);
xes(1,:) = randn(1,N);
xes(2,:) = rand(1,N) - 1/2;
xes(3,:) = cos(2*pi*7/N*n); % Prime chosen because....
xes(4,:) = cos(0.25*pi*n);
titles = ["Gaussian", "Uniformly", "cos1", "cos2"];

p = 1;
t = 1;
figure;
for i = 1:4
    x = xes(i,:);
    x_agc = AGCunity(x);
    for b = 0:4
        B = 2^b;
        xq = quant(x_agc,B);

        e = xq - x_agc;

        subplot(4,5,p);
        histogram(e);
        title(titles(t) + ", B="+B)
        p = p + 1;
    end
    t = t + 1;
end


%% lab2_4
% conclusion:
% same as in lab 2_3? i.e. signal is sampled in different points for blue signal
% but why only that weirdness in the beginning and then smooting out?
% for the other two same idea, only sampling at x points depending on 0.25/0.5
N_vals = 2.^(1:10);

figure;
for xi = 1:3
    SNR = zeros(1,length(N_vals));

    for i = 1:length(N_vals)
        N = N_vals(i);
        n = 0:N-1;
        if xi == 1
            x = cos(0.5*sqrt(3)*pi*n);
        elseif xi == 2
            x = cos(0.25*pi*n);
        elseif xi == 3
            x = cos(0.5*pi*n);
        end
        xq = quant(x,16);
        e = xq - x;
        SNR(i) = 10*log10(sum(x.^2)/sum(e.^2));
    end

    plot(N_vals, SNR); hold on;
end
ylabel("SNR (dB)")
xlabel("N");
line(N_vals, repmat(6.02*16 + 1.76, size(N_vals))); hold on;
legend("cos(0.5*sqrt(3)*pi*n)", "cos(0.25*pi*n)", "cos(0.5*pi*n)", "theoretical")

%% lab2_5
% conclusion: same as first lab
% coherent: window matters
% non-coherent: doesn't matter
N = 8192;
n = 0:N-1;
window_rect = rectwin(N)';
window_bmh = blackmanharris(N)';

x_coh = quant(cos(0.25*pi*n), 16);
x_non = quant(cos((0.25+1/N)*pi*n), 16);
X_coh_first = fft(x_coh.*window_rect, N);
X_coh_second = fft(x_coh.*window_bmh, N);
X_non_first = fft(x_non.*window_rect, N);
X_non_second = fft(x_non.*window_bmh, N);

figure;
plot(20*log10(2*abs(X_coh_first)/N)); hold on;
plot(20*log10(2*abs(X_coh_second)/N)); hold on;
plot(20*log10(2*abs(X_non_first)/N)); hold on;
plot(20*log10(2*abs(X_non_second)/N)); hold on;
legend("Coherent - rectangle", "Coherent - blackmanharris", "Non - rectangle", "Non - blackmanharris");
ylabel("Amplitude");
xlabel("F (Hz)");

%% lab2_6
% conclusions:
% noise floor lower and more "spetsigare" frequency.
% why? ... unsure
% Larger FFT => differentiate between closer frequencies dues to more freq bins
% Smaller FFT => The opposite
N = 8192;
N_2 = 16*8192;
n = 0:N-1;
n_2 = 0:N_2-1;
f_first = (0:(N-1))/N;
f_second = (0:(N_2-1))/N_2;

window_bmh_first = blackmanharris(N)';
window_bmh_second = blackmanharris(N_2)';

x_first = cos(0.5*sqrt(3)*pi*n);
x_second = cos(0.5*sqrt(3)*pi*n_2);
xq_first = quant(x_first, 16);
xq_second = quant(x_second, 16);
X_first = fft(xq_first.*window_bmh_first);
X_second = fft(xq_second.*window_bmh_second);

figure;
plot(f_first, 20*log10(abs(X_first)/max(abs(X_first)))); hold on;
plot(f_second, 20*log10(abs(X_second)/max(abs(X_second))));
legend("N = 8192", "N = 16*8192");


%% lab2_1
function y = AGCunity(x)
    y = x/max(abs(x));
end

%% lab2_2
function xq = quant(x, B)
    Q = 2^(-(B-1));
    xq = Q*(floor((x*(1-1e-10))/Q) + 1/2);

   % L = 2^B;
   % vals = Q/2*(-L+1:2:L-1);

   % xq = zeros(1, length(x));
   % for i = 1:length(x)
   %     errors = vals - x(i);
   %     [~, closest] = min(abs(errors));
   %     xq(i) = vals(closest);
   % end
end

% old
%   xq = zeros(1, length(x));
%   for i = 1:length(x)-1
%       for v = 1:length(vals)
%           e = vals(v) - x(i);
%           if -Q/2 <= abs(e) && abs(e) <= Q/2
%               xq(i) = vals(v);
%               break;
%           end
%       end
%   end
