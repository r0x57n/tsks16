clear;

N = 64;
S = [6:26, 38:57];
Q = 64;
P = 16*N;
%% about 16 to get OK results

powers = zeros(1, N);

% signals
for p = 1:P
    X = zeros(1, N);
    X(S) = qammod(randi([0 Q-1], length(S), 1)', Q);

    x = ifft(X, N);
    X_dft = fft(x);

    X_power = abs(fft(x).^2);
    X_power = X_power / max(X_power);

    powers = powers + abs(X_dft).^2;
end

powers = powers / P;
powers_db = 10*log10(powers/max(powers));

% filter
wT = 2*pi*(0:P-1)/P;
Hk = @(k) exp(1i*2*pi*k*(0:N-1)/N) / N;
expected = zeros(1, N);
for k = S
    expected = expected + abs(fft(Hk(k))).^2;
end
expected_db = 10*log10(expected/max(expected));

% plotting
f = (0:N-1)*(1/N);
figure;
hold on;
grid on;
plot(f, powers_db);
plot(f, expected_db);
legend("Simulated", "Expected")
xlabel("Normalized frequency");
ylabel("Magnitude (dB)");
