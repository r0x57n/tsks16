%clear; close all; clc;
clear; clc;

%% Part 1
L1 = 2^13;
L2 = 2^12;

M = 100;
fs = 30e6;
T = 1/fs;

F_c = 350e6;
omegac = 2*pi*F_c;

rolloff = 1/3;
S_test = 1:2:5;
Q = 64;

x1 = qammod(randi([0 Q-1], L1, 1), Q)';
x2 = qammod(randi([0 Q-1], L2, 1), Q)';

M1 = 4;
M2 = 8;
A1 = 1;
A2 = 1;
omega1 = -pi/3;
omega2 = pi/6;

c = 0.25*exp(1i*0.1*pi)*[1 zeros(1, 15) 2.4 zeros(1, 15) 1];

SNR_test = 5:20;
SINDR_x1 = zeros(1, length(SNR_test));
SINDR_x2 = zeros(1, length(SNR_test));

% S = 5 doesn't degrade to bad
% C introduces phase shifting
for S = 5
    G1 = rcosdesign(rolloff, S, M1, "sqrt")/sqrt(M1);
    G2 = rcosdesign(rolloff, S, M2, "sqrt")/sqrt(M2);
    H1 = M1*G1;
    H2 = M2*G2;

    %% Transmitter
    x1_tx = upsample(x1, M1);
    x1_tx = conv(x1_tx, H1) .* A1;
    x1_tx = x1_tx .* exp(1i*omega1*(0:length(x1_tx)-1));

    x2_tx = upsample(x2, M2);
    x2_tx = conv(x2_tx, H2) .* A2;
    x2_tx = x2_tx .* exp(1i*omega2*(0:length(x2_tx)-1));

    zeros_appended = length(x2_tx) - length(x1_tx);
    x1_tx = [x1_tx, zeros(1,zeros_appended)];
    y_tx = x1_tx + x2_tx;

    %% Zero-If
    G_lp = rcosdesign(rolloff, S, M, "sqrt")/sqrt(M);
    H_lp = M*G_lp;

        y_rx = upsample(y_tx, M);
        y_rx = conv(y_rx, H_lp);
        y_rx = y_rx .* (sqrt(2)*exp(1i*omegac/(M*fs)*(0:length(y_rx)-1)));

        y_rx = real(y_rx);
        y_rx_c = conv(y_rx, c);

    for SNR = SNR_test
        y_rx = awgn(y_rx_c, SNR, "measured");

        y_rx = y_rx .* (sqrt(2)*exp(-1i*omegac/(M*fs)*(0:length(y_rx)-1)));
        y_rx = conv(y_rx, G_lp);
        y_rx = downsample(y_rx, M);
        y_rx = y_rx(S+1:end-S);

        % Equalizer
        y_tx_ext = [y_tx, zeros(1, length(y_rx)-length(y_tx))];
        [h_eq, d_min, error_min, y_rx_eq] = eq(y_tx_ext, y_rx, M);
        y_rx_eq = y_rx_eq(d_min+1:length(y_tx)+d_min);

        %% debug
        %figure;
        %%y_rx_eq = y_rx_eq*15;
        %plot(real(y_tx), marker="O"); hold on;
        %plot(real(y_tx_ext)); hold on;
        %plot(real(y_rx)); hold on;
        %plot(real(y_rx_eq)*20); hold on;
        %plot(real(y_rx_eq_save)); hold on;
        %xlim([0 100]);
        %%ylim([-1 1])
        %legend("y_{tx}", "y_{tx_{ext}}", "y_{rx}", "y_{rx_{eq}}", "y_{rx_{eq_{save}}}");

        %% Receiver
        x1_rx = y_rx_eq .* exp(-1i*omega1*(0:length(y_rx_eq)-1));
        x1_rx = conv(x1_rx, G1);
        x1_rx = downsample(x1_rx, M1) .* 1/A1;
        x1_est = x1_rx(S+1:end-S*2);

        x2_rx = y_rx_eq .* exp(-1i*omega2*(0:length(y_rx_eq)-1));
        x2_rx = conv(x2_rx, G2);
        x2_rx = downsample(x2_rx, M2) .* 1/A2;
        x2_est = x2_rx(S+1:end-S);

        SINDR_x1(SNR) = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est - x1).^2));
        SINDR_x2(SNR) = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est - x2).^2));
    end

    %% Plotting/measuring
    figure;
    subplot(2,1,1); % scatter x1
    hold on;
    grid on;
    scatter(real(x1_est), imag(x1_est), "filled");
    scatter(real(x1), imag(x1), "filled");
    legend("x_{1est}", "x_1");

    subplot(2,1,2); % scatter x2
    hold on;
    grid on;
    scatter(real(x2_est), imag(x2_est), "filled");
    scatter(real(x2), imag(x2), "filled");
    legend("x_{2est}", "x_2");

    figure;
    hold on;
    grid on;
    plot(1:SNR_test(end), SINDR_x1);
    plot(1:SNR_test(end), SINDR_x2);
end

function [h_eq, d_min, error_min, y_rx_eq] = eq(y_tx, y_rx, N_eq)
    L = length(y_tx);
    d_min = -1;
    error_min = Inf;
    y_rx_eq = zeros(1, L);

    % toeplitz / filter
    col = [y_rx     zeros(1, N_eq)];
    row = [y_rx(1)  zeros(1, N_eq)];
    A = toeplitz(col, row);
    B = inv(A'*A)*A';

    for d = 0:N_eq
        y_tx_d = [zeros(1,d) y_tx zeros(1,N_eq-d)];
        h_eq_tmp = B*y_tx_d.';

        % eqaulized
        y_rx_eq_tmp = (A * h_eq_tmp).';

        err = sum(abs(y_rx_eq_tmp - y_tx_d).^2)/length(y_rx);

        if err < error_min
            h_eq = h_eq_tmp;
            d_min = d;
            error_min = err;
            y_rx_eq = y_rx_eq_tmp;
        end
    end
end
