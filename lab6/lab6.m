clear; clc;

%% 1)
% w0T introduces a continouse phase shift to the samples
% alpha introduces a constant phase shift of the samples
% longer signal leads to a circle, i.e. full phase rotation
% ^ leading to worse SDR

w0T = 5*pi*10^(-5);
%w0T = 0;
alpha = 0.1*pi;
%alpha = 0;

for T = 1:3
    L_set = 2.^(9:14);
    Q = 16;
    i = 1;

    if T == 1
        SDR_w0t = zeros(1, length(L_set));
        w0T = 5*pi*10^(-5);
        alpha = 0;
    elseif T == 2
        SDR_alpha = zeros(1, length(L_set));
        w0T = 0;
        alpha = 0.1*pi;
    elseif T == 3
        SDR_both = zeros(1, length(L_set));
        w0T = 5*pi*10^(-5);
        alpha = 0.1*pi;
    end

    for L = L_set
        x = qammod(randi([0 Q-1], L, 1)', Q);
        n = 1:length(x);
        y = x.*exp(1j*(w0T*n + alpha));

        %if T
        %    SDR_w0t(i) = 10*log10(sum(abs(x).^2)/sum(abs(y - x).^2));
        %elseif T == 2
        %    SDR_alpha(i) = 10*log10(sum(abs(x).^2)/sum(abs(y - x).^2));
        %elseif T == 3
        %    SDR_both(i) = 10*log10(sum(abs(x).^2)/sum(abs(y - x).^2));
        %end

        figure;
        hold on;
        scatter(real(y), imag(y), "r");
        scatter(real(x), imag(x), "b", "filled");
        grid on;
        legend("y", "x");
        title("L = " + L + ", w0T = " + w0T + ", alpha = " + alpha);


        i = i + 1;
    end
end
return;

figure;
hold on;
plot(9:14, SDR_w0t);
plot(9:14, SDR_alpha);
plot(9:14, SDR_both);
title("Q6.2) SDR");
legend("w0t", "alpha", "both");
ylabel("Mag");
xlabel("l");
return;

%% 2)
% N = 2^7 seems to work reliably for default values, any lower and it
% it depends on run to run?
% I need N = 2^8 for 30 dB SDR.
% Equalizer seems to work just as well with N_eq = 2 as with N_eq = 20, something wrong?
w0T = 5*pi*10^(-5);
alpha = 0.1*pi;

Q = 64;
L = 2^12;
SNR = 30;

tests = 1:10;
N_exps = 1:11;
SNDR = zeros(1, length(N_exps));
i = 1;
for N_exp = N_exps
    for test = tests
        SNDR_temp = zeros(1, length(N_exps));
        %N = 2^8;
        N = 2^N_exp;
        x_init = qammod(randi([0 Q-1], L, 1), Q)';

        x = [x_init(1:N), x_init(1:end)];
        y = x.*exp(1j*(w0T*(1:length(x)) + alpha));
        y = awgn(y, SNR, "measured");

        w0T_est = angle(sum(conj(y(1:N)) .* y(N+1:2*N)))/N;

        y_est = y .* exp(-1j*(w0T_est*(1:length(y))));
        [h_eq, d_min, error_min, y_eq] = eq(x, y_est, 10);
        y_comp = y_eq(d_min+1:length(x)+d_min);

        SNDR_temp(test) = 10*log10(sum(abs(x).^2)/sum(abs(y_comp - x).^2));

        % debug
        %figure;
        %plot(abs(x)); hold on;
        %plot(abs(y)); hold on;
        %plot(abs(y_eq));
        %legend("x", "y", "y_{eq}");
        %xlim([0 N*3]);
        %figure;
        %hold on;
        %scatter(real(y), imag(y), "filled");
        %scatter(real(y_est), imag(y_est), "filled");
        %scatter(real(y_comp), imag(y_comp), "filled");
        %scatter(real(x), imag(x), "filled");
        %grid on;
        %legend("y", "y_{est}", "y_{comp}", "x");
        %title("L = " + L + ", w0T = " + w0T + ", alpha = " + alpha);
        %return;
    end

    SNDR(i) = sum(SNDR_temp);

    i = i + 1;
end

figure;
plot(SNDR);
title("Q6.2) SNDR");
ylabel("Mag");
xlabel("log_2(N)");

%function [h_eq, d_min, error_min, y_rx_eq] = eq(y_tx, y_rx, N_eq)
%    L = length(y_tx);
%    d_min = -1;
%    error_min = Inf;
%    y_rx_eq = zeros(1, L);
%
%    % toeplitz / filter
%    col = [y_rx     zeros(1, N_eq)];
%    row = [y_rx(1)  zeros(1, N_eq)];
%    A = toeplitz(col, row);
%    B = inv(A'*A)*A';
%
%    for d = 0:N_eq
%        y_tx_d = [zeros(1,d) y_tx zeros(1,N_eq-d)];
%        h_eq_tmp = B*y_tx_d.';
%
%        % eqaulized
%        y_rx_eq_tmp = (A * h_eq_tmp).';
%
%        err = sum(abs(y_rx_eq_tmp - y_tx_d).^2)/length(y_rx);
%
%        if err < error_min
%            h_eq = h_eq_tmp;
%            d_min = d;
%            error_min = err;
%            y_rx_eq = y_rx_eq_tmp;
%        end
%    end
%end
