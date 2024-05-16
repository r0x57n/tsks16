clear; close all; clc;

% setup
Fs = 40e6;
Ts = 1/Fs;
delta = 1/8;
f_c1 = 10e6;
f_c2 = 30.5e6;

A_1_vals = [1, 0.01];
A_2 = 1;
A_max_vals = [0.1, 0.01];
P_vals = 0:60;

N = 8192;
x_a1 = zeros(1, N);
x_a1_phase_first = zeros(1, N);
x_a1_phase_second = zeros(1, N);
x_a2 = zeros(1, N);
t = (0:N-1)*Ts;
SIDR = zeros(1, length(P_vals));
SIDR_phase_first = zeros(1, length(P_vals));
SIDR_phase_second = zeros(1, length(P_vals));

% create the desired sequence
for k = -5:5
    f_1k = f_c1 + f_c1*k*delta;
    x_a1 = x_a1 + sin(2*pi*f_1k*t);
    x_a1_phase_first = x_a1_phase_first + sin(2*pi*f_1k*t + 0.1*(-1)^k);
    x_a1_phase_second = x_a1_phase_second + sin(2*pi*f_1k*t + 0.01*(-1)^k);
end

% save for later
y_saved_i = 1;
y_save_p = [10 20 40 60];
y_max_saved = length(y_save_p);
y_saved = zeros(y_max_saved, N);
y_saved_1 = zeros(y_max_saved, N);
y_saved_2 = zeros(y_max_saved, N);

% filter
omega_c = 0.8125*pi*Fs;
H = @(omega, P, epsilon) 1./sqrt(1 + (epsilon^2)*((omega/omega_c).^(2*P)));

legends = {};

% Plot SIDR vs P for different amplitudes
i = 1;
for A_1 = A_1_vals
    for A_max = A_max_vals
        epsilon = sqrt(10^(0.1*A_max) - 1);

        % filter the signal for multiple P:s and save the SIDR for plotting later
        for P = P_vals
            % no phase distortion
            x_a1_filtered = zeros(1,N);
            x_a2_filtered = zeros(1,N);

            % phase distortion
            x_a1_filtered_phase_first = zeros(1,N);
            x_a1_filtered_phase_second = zeros(1,N);

            for k = -5:5
                % normal
                f_1k = f_c1 + f_c1*k*delta;
                mag1 = H(2*pi*f_1k, P, epsilon);
                x_a1_filtered = x_a1_filtered + mag1*sin(2*pi*f_1k*t);

                % phase one
                x_a1_filtered_phase_first  = x_a1_filtered_phase_first + mag1*sin(2*pi*f_1k*t + 0.1*(-1)^k);
                x_a1_filtered_phase_second = x_a1_filtered_phase_second + mag1*sin(2*pi*f_1k*t + 0.01*(-1)^k);

                f_2k = f_c2 + f_c1*k*delta;
                mag2 = H(2*pi*f_2k, P, epsilon);
                x_a2_filtered = x_a2_filtered + mag2*sin(2*pi*f_2k*t);
            end

            y = A_1*x_a1_filtered + A_2*x_a2_filtered;
            y_phase_first = A_1*x_a1_filtered_phase_first + A_2*x_a2_filtered;
            y_phase_second = A_1*x_a1_filtered_phase_second + A_2*x_a2_filtered;
            y_ideal = A_1*x_a1;

            % save filter output for use later
            if (ismember(P, y_save_p))
                y_saved(y_saved_i,:) = y;
                y_saved_1(y_saved_i,:) = y_phase_first;
                y_saved_2(y_saved_i,:) = y_phase_second;
                y_saved_i = y_saved_i + 1;
            end

            num = sum(y_ideal.^2);
            den = sum((y - y_ideal).^2);
            SIDR(i) = 10*log10(num/den);
            den = sum((y_phase_first - y_ideal).^2);
            SIDR_phase_first(i) = 10*log10(num/den);
            den = sum((y_phase_second - y_ideal).^2);
            SIDR_phase_second(i) = 10*log10(num/den);

            i = i + 1;
        end

        legends{end + 1} = "A_1 = " + A_1 + ", A_{max} = " + A_max;

        % plot the SIDR vs P
        subplot(3,1,1);
        plot(P_vals, SIDR); hold on;
        grid on;
        title("Phase distortion \Delta = 0");
        legend(legends);
        ylim([-100 150]);
        ylabel("SIDR (dB)");
        xlabel("P");

        subplot(3,1,2);
        plot(P_vals, SIDR_phase_first); hold on;
        grid on;
        title("Phase distortion \Delta = 0.1");
        legend(legends);
        ylim([-100 150]);
        ylabel("SIDR (dB)");
        xlabel("P");

        subplot(3,1,3);
        plot(P_vals, SIDR_phase_second); hold on;
        grid on;
        title("Phase distortion \Delta = 0.01");
        legend(legends);
        ylim([-100 150]);
        ylabel("SIDR (dB)");
        xlabel("P");

        i = 1;
    end
end

% plot filter mag response
figure;
f = linspace(0, Fs, 2*N);
epsilon = sqrt(10^(0.1*A_max(1)) - 1);
plot(f, H(2*pi*f, 60, epsilon)); hold on;
plot(f, H(2*pi*f, 40, epsilon)); hold on;
plot(f, H(2*pi*f, 20, epsilon)); hold on;
plot(f, H(2*pi*f, 10, epsilon)); hold on;
xlim([1.5e7 4e7]);
title("Filter Magnitude Response");
legend(["P = 60", "P = 40", "P = 20", "P = 10"]);
ylabel("Mag (dB)");
xlabel("Freq (Hz)");

% plot amplitude spectrum
figure;
window = blackmanharris(N)';
f = (0:N-1)*(Fs/N);
F = (-N/2:N/2-1)*(Fs/N);
for i = 1:4
    y_w = y_saved(i,:).*window;
    Y = fft(y_w, N);
    Y_mag_db = 20*log10(abs(Y)/max(abs(Y)));
    subplot(3,1,1);
    plot(F, Y_mag_db); hold on;
    legend(num2str(y_save_p'))
    title("Phase distortion \Delta = 0");
    xlabel("Freq (Hz)");
    ylabel("Mag (dB)");

    y_w_1 = y_saved_1(i,:).*window;
    Y_1 = fft(y_w_1, N);
    Y_mag_db_1 = 20*log10(abs(Y_1)/max(abs(Y_1)));
    subplot(3,1,2);
    plot(F, Y_mag_db_1); hold on;
    legend(num2str(y_save_p'))
    title("Phase distortion \Delta = 0.1");
    xlabel("Freq (Hz)");
    ylabel("Mag (dB)");

    y_w_2 = y_saved_2(i,:).*window;
    Y_2 = fft(y_w_2, N);
    Y_mag_db_2 = 20*log10(abs(Y_2)/max(abs(Y_2)));
    subplot(3,1,3);
    plot(F, Y_mag_db_2); hold on;
    legend(num2str(y_save_p'))
    title("Phase distortion \Delta = 0.01");
    xlabel("Freq (Hz)");
    ylabel("Mag (dB)");
end




