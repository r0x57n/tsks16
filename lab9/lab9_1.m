clear;

%% vars
N = 8192;
wcT = 0.125*pi;
delta = 0.001;

%% phase noise
e1 = randn(1, N);
e2 = randn(1, N);
phi1 = zeros(1, N);
phi2 = zeros(1, N);

phi1(1) = delta*e1(1);
phi2(1) = delta*e2(1);
for n = 2:N
    phi1(n) = phi1(n-1) + delta*e1(n);
    phi2(n) = phi2(n-1) + delta*e2(n);
end

%% signal
n = 0:N-1;
x = cos(wcT*n + phi1) - 1i*sin(wcT*n + phi2);

% window
w = blackmanharris(N)';
x_w = x.*w;
X = fft(x_w, N);

%% calculate for e1 = e2
phi = zeros(1, N);
phi(1) = delta*e1(1);
for n = 2:N
    phi(n) = phi(n-1) + delta*e1(n);
end

x_equal = cos(wcT*n + phi) - 1i*sin(wcT*n + phi);
x_equal_w = x_equal.*w;
X_equal = fft(x_equal_w, N);

%% Plotting
f = linspace(0, 1, N);

figure;
subplot(2,2,1);
plot(f, 20*log10(abs(fftshift(X))));
ylim([-70 100]);
title("Amplitude x(n)");
grid on;

subplot(2,2,2);
plot(f, 20*log10(abs(fftshift(X_equal))));
ylim([-70 100]);
title("Amplitude x(n), \phi_1(n) = \phi_2(n)");
grid on;

subplot(2,2,3);
plot(f, 10*log10(abs(fftshift(X)).^2));
ylim([-70 100]);
title("Power spectrum x(n)");
grid on;

subplot(2,2,4);
plot(f, 10*log10(abs(fftshift(X_equal)).^2));
ylim([-70 100]);
title("Power spectrum x(n), \phi_1(n) = \phi_2(n)");
grid on;
