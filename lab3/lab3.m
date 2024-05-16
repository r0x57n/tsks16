% conv, downsample, fftshift, find, freqz, qammod, qamdemod, rcosdesign, scatterplot, upsample
%clear; close all; clc;
clear; close all;

%% lab3_1
% How to verify magnitude response approx unity in passband?
% Bandwidth: g1 ~ 0.3 rad/samp => 143 kHz * 2, g2 ~ 0.15 => 7 kHz * 2
% Impulse response:
rolloff1 = 1/3;
rolloff2 = 1/3;
S1 = 10; % symbol order N1 = S1*M1
S2 = 10;
M1 = 4;
M2 = 8;
L1 = 2^13;
L2 = 2^12;

g1 = rcosdesign(rolloff1, S1, M1, "sqrt")/sqrt(M1);
g2 = rcosdesign(rolloff2, S2, M2, "sqrt")/sqrt(M2);

freqz(g1, L1); hold on;
freqz(g2, L2);
legend("g1", "g2");

figure;
subplot(2,2,1);
impz(g1.*g1);
title("g1g1");
subplot(2,2,2);
impz(g2.*g2);
title("g2g2");
subplot(2,2,3);
impz(downsample(M1*g1.*g1, M1));
title("downsampled by M1");
subplot(2,2,4);
impz(downsample(M2*g2.*g2, M2));
title("downsampled by M2");

%% lab3_2 & lab3_3


fs = 30e6;
T = 1/fs;
N = 1024;

M1 = 4;
M2 = 8;

n1 = (0:L1-1)*T;
n2 = (0:L2-1)*T;
n1_up = (0:L1*M1-1)*T;
n2_up = (0:L2*M2-1)*T;

A1 = 1;
A2 = 1;
omega1 = -pi/3;
omega2 = pi/6;

S1 = 1:40;

for C = 0:2
    p = 1;

    if (C == 0)
        A1 = 1;
        A2 = 1;
    elseif (C == 1)
        A1 = 1;
        A2 = 0.1;
    elseif (C == 2)
        A1 = 0.1;
        A2 = 1;
    end

    for Q = [4, 64]
        x1 = qammod(randi([0 Q-1], L1, 1)', Q);
        x2 = qammod(randi([0 Q-1], L2, 1)', Q);

        SIDR1 = zeros(1, length(S1));
        SIDR2 = zeros(1, length(S1));
        errors1 = zeros(1, length(S1));
        errors2 = zeros(1, length(S1));

        for S = S1
            G1 = rcosdesign(rolloff1, S, M1, "sqrt")/sqrt(M1);
            G2 = rcosdesign(rolloff2, S, M2, "sqrt")/sqrt(M2);
            H1 = M1*G1;
            H2 = M2*G2;

            % Transmitter
            x1_tx = upsample(x1, M1);
            x1_tx = conv(x1_tx, H1) .* A1;
            n1_up = 1:length(x1_tx);
            x1_tx = x1_tx .* exp(1i*omega1*n1_up);

            x2_tx = upsample(x2, M2);
            x2_tx = conv(x2_tx, H2) .* A2;
            n2_up = 1:length(x2_tx);
            x2_tx = x2_tx .* exp(1i*omega2*n2_up);

            zeros_appended = length(x2_tx) - length(x1_tx);
            x1_tx = [x1_tx, zeros(1,zeros_appended)];
            y_tx = x1_tx + x2_tx;

            % Receiver
            n1_up = 1:length(y_tx);
            x1_rx = y_tx .* exp(-1i*omega1*n1_up);
            x1_rx = conv(x1_rx, G1);
            x1_rx = downsample(x1_rx, M1) .* 1/A1;
            x1_est = x1_rx(S+1:end-S*2);

            n2_up = 1:length(y_tx);
            x2_rx = y_tx .* exp(-1i*omega2*n2_up);
            x2_rx = conv(x2_rx, G2);
            x2_rx = downsample(x2_rx, M2) .* 1/A2;
            x2_est = x2_rx(S+1:end-S);

            SIDR1(S) = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est - x1).^2));
            SIDR2(S) = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est - x2).^2));
            errors1(S) = length(find(qamdemod(x1,Q) - qamdemod(x1_est,Q)));
            errors2(S) = length(find(qamdemod(x2,Q) - qamdemod(x2_est,Q)));

            %scatterplot(x1)
            %title("x_1(n_1)")
            %scatterplot(x1_est)
            %title("x_{1est}")

            %scatterplot(x2)
            %title("x_2(n_2)")
            %scatterplot(x2_est)
            %title("x_{2est}")
            %return;
        end

        if (A1 == 1 && A2 == 1)
            figure(4);
        elseif (A1 == 1 && A2 == 0.1)
            figure(6);
        elseif (A1 == 0.1 && A2 == 1)
            figure(8);
        end
        subplot(2,1,p);
        plot(S1, SIDR1); hold on;
        plot(S1, SIDR2);
        xlabel("S");
        ylabel("SIDR (dB)");
        title("A_1 = " + A1 + ", A_2 = " + A2 + ", Q = " + Q);
        legend("x1", "x2")

        if (A1 == 1 && A2 == 1)
            figure(5);
        elseif (A1 == 1 && A2 == 0.1)
            figure(7);
        elseif (A1 == 0.1 && A2 == 1)
            figure(9);
        end
        subplot(2,1,p);
        plot(S1, errors1); hold on;
        plot(S1, errors2)
        xlabel("S");
        ylabel("Errors");
        title("A_1 = " + A1 + ", A_2 = " + A2 + ", Q = " + Q);
        legend("x1", "x2");

        p = p + 1;
    end
end

%% lab3_4
A1 = 1;
A2 = 0.01;

p = 1;
for Q = [4, 64]
    x1 = qammod(randi([0 Q-1], L1, 1)', Q);
    x2 = qammod(randi([0 Q-1], L2, 1)', Q);

    SIDR1 = zeros(1, length(S1));
    SIDR2 = zeros(1, length(S1));
    errors1 = zeros(1, length(S1));
    errors2 = zeros(1, length(S1));

    for S = 20
        G1 = rcosdesign(rolloff1, S, M1, "sqrt")/sqrt(M1);
        G2 = rcosdesign(rolloff2, S, M2, "sqrt")/sqrt(M2);
        H1 = M1*G1;
        H2 = M2*G2;

        % Transmitter
        x1_tx = upsample(x1, M1);
        x1_tx = conv(x1_tx, H1) .* A1;
        n1_up = 1:length(x1_tx);
        x1_tx = x1_tx .* exp(1i*omega1*n1_up);

        x2_tx = upsample(x2, M2);
        x2_tx = conv(x2_tx, H2) .* A2;
        n2_up = 1:length(x2_tx);
        x2_tx = x2_tx .* exp(1i*omega2*n2_up);

        zeros_appended = length(x2_tx) - length(x1_tx);
        x1_tx = [x1_tx, zeros(1,zeros_appended)];
        y_tx = x1_tx + x2_tx;

        % Receiver
        n1_up = 1:length(y_tx);
        x1_rx = y_tx .* exp(-1i*omega1*n1_up);
        x1_rx = conv(x1_rx, G1);
        x1_rx = downsample(x1_rx, M1) .* 1/A1;
        x1_est = x1_rx(S+1:end-S*2);

        n2_up = 1:length(y_tx);
        x2_rx = y_tx .* exp(-1i*omega2*n2_up);
        x2_rx = conv(x2_rx, G2);
        x2_rx = downsample(x2_rx, M2) .* 1/A2;
        x2_est = x2_rx(S+1:end-S);

        SIDR1(S) = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est - x1).^2));
        SIDR2(S) = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est - x2).^2));
        errors1(S) = length(find(qamdemod(x1,Q) - qamdemod(x1_est,Q)));
        errors2(S) = length(find(qamdemod(x2,Q) - qamdemod(x2_est,Q)));
    end

    figure(9);
    subplot(2,1,p);
    plot(S1, errors1); hold on;
    plot(S1, errors2)
    xlabel("S");
    ylabel("Errors");
    title("A_1 = " + A1 + ", A_2 = " + A2 + ", Q = " + Q);
    legend("x1", "x2");

    p = p + 1;
end
