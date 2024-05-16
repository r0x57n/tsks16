clear; clc;

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
SNR = 20;
S = 5;

x1 = qammod(randi([0 Q-1], L1, 1)', Q);
x2 = qammod(randi([0 Q-1], L2, 1)', Q);

M1 = 4;
M2 = 8;
A1 = 1;
A2 = 1;
omega1 = -pi/3;
omega2 = pi/6;

w0T = 5*pi*10^(-5);
alpha = 0.1*pi;
%N = 128;
%N = 256;
N = 512;

c = 0.25*exp(1i*0.1*pi)*[1 zeros(1, 15) 2.4 zeros(1, 15) 1];

G1 = rcosdesign(rolloff, S, M1, "sqrt")/sqrt(M1);
G2 = rcosdesign(rolloff, S, M2, "sqrt")/sqrt(M2);
H1 = M1*G1;
H2 = M2*G2;

%%
%% Transmitter
%%
x1_tx = upsample(x1, M1);
x1_tx = conv(x1_tx, H1) .* A1;
x1_tx = x1_tx .* exp(1i*omega1*(0:length(x1_tx)-1));

x2_tx = upsample(x2, M2);
x2_tx = conv(x2_tx, H2) .* A2;
x2_tx = x2_tx .* exp(1i*omega2*(0:length(x2_tx)-1));

zeros_appended = length(x2_tx) - length(x1_tx);
x1_tx = [x1_tx, zeros(1,zeros_appended)];

y_tx = x1_tx + x2_tx;

% extend (for CFO/PO)
y_tx = [y_tx(1:N), y_tx(1:end)];

%%
%% Zero-If
%%
G_lp = rcosdesign(rolloff, S, M, "sqrt")/sqrt(M);
H_lp = M*G_lp;

% send through channel
y_rx = upsample(y_tx, M);
y_rx = conv(y_rx, H_lp);
y_rx = y_rx .* (sqrt(2)*exp(1i*omegac/(M*fs)*(0:length(y_rx)-1)));

y_rx = real(y_rx);
y_rx = conv(y_rx, c);
y_rx = awgn(y_rx, SNR, "measured");

% introduce CFO/PO
y_rx = y_rx .* (sqrt(2)*exp(-1i*((omegac - w0T)/(M*fs))*(0:length(y_rx)-1) - alpha));
y_rx = conv(y_rx, G_lp);
y_rx = downsample(y_rx, M);

y_rx = y_rx(S+1:end-S); % remove filter delay

%%
%% Compensation
%%
w0T_est = angle(sum(conj(y_rx(1:N)) .* y_rx(N+1:2*N)))/N;
%y_rx = y_rx.*exp(-1j*(w0T_est*(1:length(y_rx))));
%% should be N+1?
y_rx = y_rx(N+2:end); % remove extension

y_tx_ext = [y_tx(N+1:end) zeros(1, length(y_rx)-length(y_tx))];
[h_eq, d_min, error_min, y_eq] = eq(y_tx_ext, y_rx, M);
y_comp = y_eq(d_min+1:length(y_tx(N+1:end))+d_min);

%%
%% Receiver
%%
x1_rx = y_comp .* exp(-1i*omega1*(0:length(y_comp)-1));
x1_rx = conv(x1_rx, G1);
x1_rx = downsample(x1_rx, M1) .* 1/A1;
x1_est = x1_rx(S+1:end-S*2);

x2_rx = y_comp .* exp(-1i*omega2*(0:length(y_comp)-1));
x2_rx = conv(x2_rx, G2);
x2_rx = downsample(x2_rx, M2) .* 1/A2;
x2_est = x2_rx(S+1:end-S);

%%
%% Plotting/measuring
%%

figure;
subplot(2,1,1);
hold on;
grid on;
scatter(real(x1_est), imag(x1_est), "filled");
scatter(real(x1), imag(x1), "filled");
legend("x_{1est}", "x_1");
subplot(2,1,2);
hold on;
grid on;
scatter(real(x2_est), imag(x2_est), "filled");
scatter(real(x2), imag(x2), "filled");
legend("x_{2est}", "x_2");

x1_SNDR = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est - x1).^2));
x2_SNDR = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est - x2).^2));

x1 = qamdemod(x1, Q);
x2 = qamdemod(x2, Q);
x1_est = qamdemod(x1_est, Q);
x2_est = qamdemod(x2_est, Q);

x1_errors = sum(x1_est - x1);
x2_errors = sum(x2_est - x2);

fprintf("x1:\n")
fprintf("symbol errors: %f\n", x1_SNDR);
fprintf("SNDR: %f\n", x1_errors);

fprintf("x2:\n")
fprintf("symbol errors: %f\n", x2_SNDR);
fprintf("SNDR: %f\n\n", x2_errors);

fprintf("N = 128 is too low, N = 512 works well, N = 256 also works okayish\n");
fprintf("N = 512 gives better SDR (almost 30), not as told by first parts...\n");
fprintf("d_min is usually around 30-70 so 100 is probably good\n");

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
