% Parameters
Z0 = 50; % Characteristic impedance [Ohm]
ZL = 704; % Load impedance [Ohm]
f = 2242e6; % Frequency [Hz]
c = 3e8; % Speed of light [m/s]
lambda = c / f; % Wavelength [m]

% Reflection Coefficient at the Load
gamma_L = (ZL - Z0) / (ZL + Z0);

% Phase constant
beta = 2 * pi / lambda;

% Distance and stub length (optimized)
d = lambda / 4; % Stub distance from the load
l = lambda / 4; % Stub length (minimum length)

% Frequency range for analysis
f_min = 0.1 * f;
f_max = 1.9 * f;
frequencies = linspace(f_min, f_max, 500);

% Initialize matrices for calculations
ABCD_params = zeros(2, 2, length(frequencies));
S_params = zeros(2, 2, length(frequencies));
Gamma_left = zeros(1, length(frequencies));
Z_input_left = zeros(1, length(frequencies));
Z_input_right = zeros(1, length(frequencies));
Voltage_10mm = zeros(1, length(frequencies));
Current_10mm = zeros(1, length(frequencies));
Power_10mm = zeros(1, length(frequencies));
Reflected_power = zeros(1, length(frequencies));
Transmitted_power = zeros(1, length(frequencies));
Current_left = zeros(1, length(frequencies));
Current_right = zeros(1, length(frequencies));

% Incident Voltage Wave
V_in = 5; % Amplitude [V]

% Main loop for frequency-dependent calculations
for idx = 1:length(frequencies)
    omega = 2 * pi * frequencies(idx);
    lambda_f = c / frequencies(idx);
    beta_f = 2 * pi / lambda_f;

    % Transmission line equations for ABCD matrix
    T_stub = [1 0; 0 1]; % Placeholder for stub calculations
    T_line = [cos(beta_f * d) 1j * Z0 * sin(beta_f * d); ...
              1j * sin(beta_f * d) / Z0 cos(beta_f * d)];
    
    % ABCD matrix calculation
    ABCD_params(:, :, idx) = T_stub * T_line;

    % Updated S11 and S21 formulas
    A = ABCD_params(1, 1, idx);
    B = ABCD_params(1, 2, idx);
    C = ABCD_params(2, 1, idx);
    D = ABCD_params(2, 2, idx);
    S_params(1, 1, idx) = (A + B - C - D) / (A + B + C + D);
    S_params(1, 2, idx) = 2 / (A + B + C + D);

    % Reflection coefficient immediately on the left of the stub
    Gamma_left(idx) = S_params(1, 1, idx);

    % Input impedance immediately on the left of the stub
    Z_input_left(idx) = Z0 * (1 + Gamma_left(idx)) / (1 - Gamma_left(idx));

    % Input impedance immediately on the right of the stub
    Z_input_right(idx) = ZL;

    % Voltage and current at 10mm from the load
    distance_10mm = 0.01; % Distance from the load in meters
    voltage_phase = exp(-1j * beta_f * distance_10mm);
    Voltage_10mm(idx) = abs(V_in * voltage_phase);
    Current_10mm(idx) = Voltage_10mm(idx) / Z_input_right(idx);

    % Active power at 10mm
    Power_10mm(idx) = 0.5 * abs(Voltage_10mm(idx))^2 / real(Z_input_right(idx));

    % Reflected and transmitted power
    Reflected_power(idx) = abs(S_params(1, 1, idx))^2;
    Transmitted_power(idx) = 1 - Reflected_power(idx);

    % Stub currents
    Voltage_left = V_in * (1 + Gamma_left(idx));
    Current_left(idx) = Voltage_left / Z_input_left(idx);
    Voltage_right = V_in * (1 - Gamma_left(idx));
    Current_right(idx) = Voltage_right / Z_input_right(idx);
end

% Plot for Question 1: Scattering Coefficients
figure;
subplot(2, 1, 1);
plot(frequencies, abs(squeeze(S_params(1, 1, :))), frequencies, abs(squeeze(S_params(1, 2, :))));
title('1. Scattering Coefficients (Modulus)');
xlabel('Frequency (Hz)');
ylabel('|S|');
legend('S11', 'S21');

subplot(2, 1, 2);
plot(frequencies, angle(squeeze(S_params(1, 1, :))), frequencies, angle(squeeze(S_params(1, 2, :))));
title('1. Scattering Coefficients (Phase)');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');

% Question 2: Reflection Coefficient
% figure;
% plot(frequencies, abs(Gamma_left));
% title('2. Reflection Coefficient');
% xlabel('Frequency (Hz)');
% ylabel('|Gamma|');

% Modulus Plot
subplot(2, 1, 1);
plot(frequencies, abs(Gamma_left)); 
title('2. Modulus of Reflection Coefficient (\Gamma) Immediately Left of the Stub');
xlabel('Frequency (Hz)');
ylabel('|Gamma|');

% Phase Plot
subplot(2, 1, 2);
plot(frequencies, angle(Gamma_left)); 
title('2. Phase of Reflection Coefficient (\Gamma) Immediately Left of the Stub');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');

% Question 3: Input Impedance Left
figure;
subplot(2, 1, 1);
plot(frequencies, real(Z_input_left));
title('3. Input Impedance (Real Part, Left of Stub)');
xlabel('Frequency (Hz)');
ylabel('Real(Z) [Ohm]');

subplot(2, 1, 2);
plot(frequencies, imag(Z_input_left));
title('3. Input Impedance (Imaginary Part, Left of Stub)');
xlabel('Frequency (Hz)');
ylabel('Imag(Z) [Ohm]');

% Question 4: Input Impedance Right
figure;
subplot(2, 1, 1);
plot(frequencies, real(Z_input_right));
title('4. Input Impedance (Real Part, Right of Stub)');
xlabel('Frequency (Hz)');
ylabel('Real(Z) [Ohm]');

subplot(2, 1, 2);
plot(frequencies, imag(Z_input_right));
title('4. Input Impedance (Imaginary Part, Right of Stub)');
xlabel('Frequency (Hz)');
ylabel('Imag(Z) [Ohm]');

% Question 5: Stub Left Side Current
figure;
subplot(2, 1, 1);
plot(frequencies, abs(Current_left));
title('5. Stub Left Side Current (Modulus)');
xlabel('Frequency (Hz)');
ylabel('Current (A)');

subplot(2, 1, 2);
plot(frequencies, angle(Current_left));
title('5. Stub Left Side Current (Phase)');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');

% Question 6: Stub Right Side Current
figure;
subplot(2, 1, 1);
plot(frequencies, abs(Current_right));
title('6. Stub Right Side Current (Modulus)');
xlabel('Frequency (Hz)');
ylabel('Current (A)');

subplot(2, 1, 2);
plot(frequencies, angle(Current_right));
title('6. Stub Right Side Current (Phase)');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');

% Question 7: Reflected and Transmitted Power
figure;
subplot(2, 1, 1);
plot(frequencies, Reflected_power);
title('7. Reflected Power');
xlabel('Frequency (Hz)');
ylabel('Power Ratio');

subplot(2, 1, 2);
plot(frequencies, Transmitted_power);
title('7. Transmitted Power');
xlabel('Frequency (Hz)');
ylabel('Power Ratio');

% Question 8: Modified Stub Length Analysis
l_prime = l + lambda; % New stub length
Gamma_left_prime = zeros(1, length(frequencies));
Z_input_left_prime = zeros(1, length(frequencies));

for idx = 1:length(frequencies)
    omega = 2 * pi * frequencies(idx);
    lambda_f = c / frequencies(idx);
    beta_f = 2 * pi / lambda_f;

    % ABCD Matrix with new stub length
    T_stub_prime = [cos(beta_f * l_prime) 1j * Z0 * sin(beta_f * l_prime); ...
                   1j * sin(beta_f * l_prime) / Z0 cos(beta_f * l_prime)];
    T_line = [cos(beta_f * d) 1j * Z0 * sin(beta_f * d); ...
              1j * sin(beta_f * d) / Z0 cos(beta_f * d)];
    ABCD_prime = T_stub_prime * T_line;

    % Recalculate S11 for modified stub length
    A_prime = ABCD_prime(1, 1);
    B_prime = ABCD_prime(1, 2);
    C_prime = ABCD_prime(2, 1);
    D_prime = ABCD_prime(2, 2);
    S_params_prime_11 = (A_prime + B_prime - C_prime - D_prime) / (A_prime + B_prime + C_prime + D_prime);
    Gamma_left_prime(idx) = S_params_prime_11;

    % Recalculate input impedance left of stub for modified length
    Z_input_left_prime(idx) = Z0 * (1 + Gamma_left_prime(idx)) / (1 - Gamma_left_prime(idx));
end

% Plot for Modified Stub Length
figure;
subplot(2, 1, 1);
plot(frequencies, real(Z_input_left_prime));
title('8. Real Part of Input Impedance (Modified Stub Length)');
xlabel('Frequency (Hz)');
ylabel('Real(Z) [Ohm]');

subplot(2, 1, 2);
plot(frequencies, imag(Z_input_left_prime));
title('8. Imaginary Part of Input Impedance (Modified Stub Length)');
xlabel('Frequency (Hz)');
ylabel('Imag(Z) [Ohm]');

% Question 9: Voltage at 10mm
figure;
plot(frequencies, Voltage_10mm);
title('9. Voltage at 10mm from the Load');
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');

% Question 10: Current at 10mm
figure;
plot(frequencies, Current_10mm);
title('10. Current at 10mm from the Load');
xlabel('Frequency (Hz)');
ylabel('Current (A)');

% Question 11: Power at 10mm
figure;
plot(frequencies, Power_10mm);
title('11. Active Power at 10mm from the Load');
xlabel('Frequency (Hz)');
ylabel('Power (W)');
