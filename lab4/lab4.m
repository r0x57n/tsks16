clear; close all; clc;

L = 2^12;
N_eq = 20;
c = 0.25*exp(1i*0.1*pi)*[1 2.4 1];

%% Case 1
y_tx_1 = randn(1, L);

% rx
y_rx_1 = conv(y_tx_1, c);
y_rx_1 = awgn(y_rx_1, 30, "measured");

% calc/plot
[h_eq, d_min, error_min] = eq(y_tx_1, y_rx_1, N_eq);
fprintf("Case 1:\nd_min: %d\nerror_min: %d\n\n", d_min, error_min);
figure;
freqz(h_eq);
title("Case 1");

%% Case 2
rolloff = 0.25;
span = 6;
sps = 2;
rcosFilter = rcosdesign(rolloff, span, sps);

y_tx_2 = filter(rcosFilter, 1, y_tx_1);

% rx
y_rx_2 = conv(y_tx_2, c);
y_rx_2 = awgn(y_rx_2, 30, "measured");

% calc/plot
[h_eq, d_min, error_min] = eq(y_tx_2, y_rx_2, N_eq);
fprintf("Case 2:\nd_min: %d\nerror_min: %d\n\n", d_min, error_min);
figure;
freqz(h_eq);
title("Case 2");

%% part 3
SNDR = zeros(1, 50);
for N = 1:50
    [~, ~, error_min] = eq(y_tx_1, y_rx_1, N);
    SNDR(N) = 10*log10(sum(abs(y_tx_1).^2)/error_min);
end
figure;
plot(SNDR);
title("Case 1");

for N = 1:50
    [~, ~, error_min] = eq(y_tx_2, y_rx_2, N);
    SNDR(N) = 10*log10(sum(abs(y_tx_2).^2)/error_min);
end
figure;
plot(SNDR);
title("Case 2");

%% function
function [h_eq, d_min, error_min] = eq(y_tx, y_rx, N_eq)
    L = length(y_tx);
    d_min = -1;
    error_min = Inf;

    for d = 0:N_eq
        y_tx_d = [zeros(1,2) y_tx(1:end)];
        y_tx_d = [zeros(1,d) y_tx_d zeros(1,N_eq-d)];

        % toeplitz / filter
        col = [y_rx     zeros(1, N_eq)];
        row = [y_rx(1)  zeros(1, N_eq)];
        A = toeplitz(col, row);
        B = inv(A'*A)*A';
        h_eq_tmp = B*y_tx_d';

        % eqaulized
        y_rx_eq = conv(y_rx, h_eq_tmp);

        err = sum(abs(y_rx_eq - y_tx_d).^2);

        if err < error_min
            h_eq = h_eq_tmp;
            d_min = d;
            error_min = err;
        end
    end
end
