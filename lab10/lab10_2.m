clear;

N_values = [32, 1024, 8192];
Q = 4;
tests = 10;

figure;
set(gcf, 'Position', [50, 1000, 1000, 700])
splot = 1;

for N = N_values
    K_values = 1:N;
    PAPR_values = zeros(1, length(K_values));

    for K = K_values
        papr_sum = 0;

        for i = 1:tests
            X = zeros(1, N);
            S = randperm(N, K);
            X(S) = qammod(randi([0 Q-1], 1, K), Q, 'UnitAveragePower', true);

            X = X / sqrt(mean(abs(X).^2));

            x = ifft(X) * sqrt(N);

            power = abs(x).^2;
            PAPR = max(power) / mean(power);
            papr_sum = papr_sum + PAPR;
        end

        PAPR_values(K) = papr_sum / tests;
    end

    if Q == 64
        theo = 10 * log10(7 * K_values / 3);
    else
        theo = 10 * log10(K_values);
    end

    subplot(3, 1, splot);
    hold on;
    grid on;
    plot(K_values, 10*log10(PAPR_values), "b");
    plot(K_values, theo, "r--");
    legend("Simulated", "Theoretical");
    xlabel("Active Carriers (K)");
    ylabel("PAPR (dB)");
    title("PAPR for N = " + N + ", Q = " + Q + ", tests = " + tests);

    splot = splot + 1;
end
