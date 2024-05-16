clear;

%%
%% Variables
%%

%% Signals
Q = 16;
L1 = 2^13;
L2 = 2^12;
fs = 30e6;
T = 1/fs;
x1 = qammod(randi([0 Q-1], L1, 1)', Q);
x2 = qammod(randi([0 Q-1], L2, 1)', Q);

%% Misc
c = 0.25*[1 zeros(1, 15) 2.4 zeros(1, 15) 1];
rolloff = 1/3;
M = 100;

%% Transmitter/Receiver
M1 = 4;
M2 = 8;
%A1 = 0.05;
A1 = 0;
A2 = 1;
S1 = 30;
S2 = 30;
omega1 = -pi/3;
omega2 = pi/6;
G1 = rcosdesign(rolloff, S1, M1, "sqrt")/sqrt(M1);
G2 = rcosdesign(rolloff, S2, M2, "sqrt")/sqrt(M2);
H1 = M1*G1;
H2 = M2*G2;

%% Zero-If
SNR = 50;
S = 20;
F_c = 350e6;
wc = 2*pi*F_c;

%% CFO/PO
w0T = 5*pi*10^(-5);
alpha = 0.1*pi;
N = 512;

%% IQM & SFM
g1 = 0.99;
g2 = 1.02;
phi1 = 0.5;
phi2 = 0.55;
%g1 = 0;
%g2 = 0;
%phi1 = 0;
%phi2 = 0;
r = 3*10^(-5);

%% Equalizer
N_eq = 6;

%% Nonlinearities
F_p1 = -10e6;
F_p2 = 10e6;
F_s1 = -15e6;
F_s2 = 15e6;
delta_p = 0.001;
delta_s = 0.001;

f_p1 = (F_p1 + fs/2) / fs;
f_p2 = (F_p2 + fs/2) / fs;
f_s1 = (F_s1 + fs/2) / fs;
f_s2 = (F_s2 + fs/2) / fs;

[BP_N, fo, ao, w] = firpmord([f_s1 f_p1 f_p2 f_s2], [0 1 0], [delta_s delta_p delta_s]);
BP = firpm(BP_N, fo, ao, w);

%% Errors / Fixes

% Errors
CHANNEL     = 1;
AWGN        = 1;
IQM         = 1;
SFM         = 0; % only part of Lab 7 last task?
NONLIN      = 1;

% Fixes
NONLIN_FIX  = 1;
NONLIN_BP   = 1;
IQM_FIX     = 1;
CFO_PO      = 1;
CFO_FIX     = 1;
EQ          = 1;

%%
%% Transmitter
%%
x1_tx = upsample(x1, M1);
x1_tx = conv(x1_tx, H1) .* A1;
x1_m = 0:length(x1_tx) - 1;
x1_tx = x1_tx .* exp(1i*omega1*x1_m);

x2_tx = upsample(x2, M2);
x2_tx = conv(x2_tx, H2) .* A2;
x2_m = 0:length(x2_tx) - 1;
x2_tx = x2_tx .* exp(1i*omega2*x2_m);

zeros_appended = length(x2_tx) - length(x1_tx);
x1_tx = [x1_tx, zeros(1,zeros_appended)];

y_tx = x1_tx + x2_tx;

% extend (for CFO/PO)
if CFO_PO
    y_tx = [y_tx(1:N), y_tx];
end

%%
%% Zero-If
%%
G_lp = rcosdesign(rolloff, S, M, "sqrt")/sqrt(M);
H_lp = M*G_lp;

% Transmit
y_rx = upsample(y_tx, M);
y_rx = conv(y_rx, H_lp);
y_rx_m = 0:length(y_rx) - 1;
y_rx = y_rx .* (sqrt(2)*exp(1i*wc/(M*fs)*y_rx_m));

y_rx = real(y_rx);
if CHANNEL
    y_rx = conv(y_rx, c);
end

if AWGN
    y_rx = awgn(y_rx, SNR, "measured");
end

if NONLIN
    a0 = 0.01;
    a1 = 1;
    a2 = -0.2/max(abs(y_rx));
    a3 = 0.15/max(abs(y_rx).^2);
    y_rx = a0 + a1.*y_rx + a2.*(y_rx).^2 + a3.*(y_rx).^3;
end

if NONLIN_BP
    y_rx = conv(y_rx, BP);
end

% Demodulator

% be bad to phase (IQM/CFO/PO)
y_rx_m = 0:length(y_rx) - 1;
Th = 1/(M*fs);
phase = exp(-1i*wc*Th*y_rx_m);
if IQM && CFO_PO
    phase = g1*cos(((wc - w0T)*Th)*y_rx_m + phi1 - alpha) - 1i*g2*sin(((wc - w0T)*Th)*y_rx_m + phi2 - alpha);
elseif IQM && ~CFO_PO
    phase = g1*cos((wc*Th)*y_rx_m + phi1) - 1i*g2*sin((wc*Th)*y_rx_m + phi2);
elseif ~IQM && CFO_PO
    fprintf("apply2\n");
    phase = exp(-1i*(wc - w0T)*Th*y_rx_m - alpha);
end

y_rx = y_rx.*(sqrt(2)*phase);

% Receive
y_rx = conv(y_rx, G_lp);
y_rx = downsample(y_rx, M);

% remove filter delay
%% added this -1 because of the channel, should it be there?
if CHANNEL
    % if no channel, let equalizer fix mess!
    % actually, always let equalizer fix mess, and if there's
    % no channel, no need to fix anything
    y_rx = y_rx(S+1:end-S-1);
else
    y_rx = y_rx(S+1:end-S);
end

% introduce SFM
if SFM
    y_rx = Interpolation_Farrow(y_rx, 1+r);
end

%%
%% Compensation
%%

% order:
%    estimate IQM
%    compensate IQM
%    estimate CFO
%    compensate CFO
%    remove the double samples
%    then equalize

y_comp = y_rx;

if NONLIN_FIX
    y_comp = y_comp - a3.*(abs(y_comp).^2).*y_comp;
end

if IQM
    % IQM estimation
    yr = real(y_comp);
    yi = imag(y_comp);
    Pr = mean(yr.^2);
    Pi = mean(yi.^2);
    g_est = sqrt(Pi/Pr);
    yp = -yr.*yi;
    Pyp = mean(yp);
    phi_est = Pyp/Pr;

    % IQM Compensation
    yr = yr*g_est;
    yi = yi + yr*phi_est;

    if IQM_FIX
        y_comp = yr + 1i*yi;
    end
end

% CFO/PO
if CFO_PO
    % estimate
    w0T_est = angle(sum(conj(y_comp(1:N)) .* y_comp(N+1:2*N)))/N;

    % compensate
    if CFO_FIX
        y_comp = y_comp.*exp(-1j*(w0T_est*(1:length(y_comp))));
    end

    % remove extension
    y_comp = y_comp(N+1:end);
end

% Equalizing
if EQ
    %% With CFO/PO
    if CFO_PO
        y_tx_ext = [y_tx(N+1:end) zeros(1, length(y_comp)-length(y_tx))];
        y_eq_end = length(y_tx(N+1:end));
    %% Without CFO/PO
    elseif ~CFO_PO
        y_eq_end = length(y_tx);
        y_tx_ext = [y_tx zeros(1, length(y_comp)-length(y_tx))];
    end

    [h_eq, d_min, error_min, y_eq] = eq(y_tx_ext, y_comp, N_eq);
    y_comp = y_eq(d_min+1 : y_eq_end+d_min);
end

%%
%% Receiver
%%
m = 0:length(y_comp) - 1;

x1_rx = y_comp .* exp(-1i*omega1*m);
x1_rx = conv(x1_rx, G1);
x1_rx = downsample(x1_rx, M1) .* 1/A1;
x1_est = x1_rx(S1+1:end-S1*2);

x2_rx = y_comp .* exp(-1i*omega2*m);
x2_rx = conv(x2_rx, G2);
x2_rx = downsample(x2_rx, M2) .* 1/A2;
x2_est = x2_rx(S2+1:end-S2);

%%
%% Plotting/measuring
%%

figure;
set(gcf,'position',[50, 1000, 1000, 700])
subplot(3,2,1); % scatter x1
xlim([-5 5]);
ylim([-5 5]);
hold on;
grid on;
scatter(real(x1_est), imag(x1_est), "filled");
scatter(real(x1), imag(x1), "filled");
legend("x_{1est}", "x_1");

subplot(3,2,2); % scatter x2
xlim([-5 5]);
ylim([-5 5]);
hold on;
grid on;
scatter(real(x2_est), imag(x2_est), "filled");
scatter(real(x2), imag(x2), "filled");
legend("x_{2est}", "x_2");

%% spectrum stuff
%subplot(3,2,3); % spectrum x1
%hold on;
%grid on;
%pspectrum(x1_est);
%pspectrum(x1);
%legend("x1_{est}", "x1");

%subplot(3,2,4); % spectrum x2
%hold on;
%grid on;
%pspectrum(x2_est);
%pspectrum(x2);
%legend("x2_{est}", "x2");

%% signals
n_comp = length(y_comp);
n_y_tx = length(y_tx);
n_y_rx = length(y_rx);

f_comp = (-n_comp/2:n_comp/2-1)*(fs/n_comp);
f_y_tx = (-n_y_tx/2:n_y_tx/2-1)*(fs/n_y_tx);
f_y_rx = (-n_y_rx/2:n_y_rx/2-1)*(fs/n_y_rx);

Y_comp = fftshift(fft(y_comp));
Y_y_tx = fftshift(fft(y_tx));
Y_y_rx = fftshift(fft(y_rx));

P_comp = abs(Y_comp).^2 / n_comp;
P_y_tx = abs(Y_y_tx).^2 / n_y_tx;
P_y_rx = abs(Y_y_rx).^2 / n_y_rx;

subplot(3,2,3);
ylim([-150 50]);
plot(f_y_tx, 10*log10(P_y_tx));
grid on;
legend("y_{tx}");

subplot(3,2,4);
ylim([-150 50]);
hold on;
plot(f_y_rx, 10*log10(P_y_rx));
grid on;
legend("y_{rx}");

subplot(3,2,5);
plot(sin(x1));
legend(":)");
grid on;

subplot(3,2,6);
ylim([-150 50]);
hold on;
plot(f_comp, 10*log10(P_comp));
grid on;
legend("y_{comp}");

%% values
%if NONLIN_BP
%    x1_SNDR = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est(1:end-1) - x1).^2));
%    x2_SNDR = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est(1:end-1) - x2).^2));
%else
    x1_SNDR = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est(1:end) - x1).^2));
    x2_SNDR = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est(1:end) - x2).^2));
%end

x1 = qamdemod(x1, Q);
x2 = qamdemod(x2, Q);
x1_est = qamdemod(x1_est, Q);
x2_est = qamdemod(x2_est, Q);

%if NONLIN_BP
%    x1_errors = sum(x1_est(1:end-1) - x1);
%    x2_errors = sum(x2_est(1:end-1) - x2);
%else
    x1_errors = sum(x1_est - x1);
    x2_errors = sum(x2_est - x2);
%end

fprintf("x1:\n")
fprintf("symbol errors: %f\n", x1_errors);
fprintf("SNDR: %f\n", x1_SNDR);

fprintf("x2:\n")
fprintf("symbol errors: %f\n", x2_errors);
fprintf("SNDR: %f\n\n", x2_SNDR);

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
