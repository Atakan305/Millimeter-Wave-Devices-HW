% Constants
c = 3e8; % Speed of light in vacuum (m/s)
a = 0.0740; % Waveguide width (m)
b = 0.0296; % Waveguide height (m)
epsilon_r = 3.7; % Relative permittivity
epsilon_0 = 8.854e-12; % Vacuum permittivity coefficient (F/m)
mu_0 = 4 * pi * 1e-7; % Vacuum permeability coefficient (H/m)

% Mode numbers
m = 3; 
n = 4;

% Cutoff frequency
f_c = (c / 2) * sqrt((m / a)^2 + (n / b)^2) / sqrt(epsilon_r);

% Frequency range
f = linspace(1.3 * f_c, 1.9 * f_c, 500); % Interval from 1.3fc to 1.9fc
omega = 2 * pi * f;

% Cutoff wavenumber
k_c = sqrt((m * pi / a)^2 + (n * pi / b)^2);

% Initialize arrays
beta = zeros(size(f));
v_g = zeros(size(f));
v_f = zeros(size(f));
Z_TM = zeros(size(f));
lambda_g = zeros(size(f));

% Loop through frequencies
for i = 1:length(f)
    k = omega(i) / c; % Wavenumber
    beta_squared = k^2 - k_c^2;
    if beta_squared > 0
        beta(i) = sqrt(beta_squared);
        v_g(i) = c * sqrt(1 - (f_c / f(i))^2);
        v_f(i) = c / sqrt(1 - (f_c / f(i))^2);
        Z_TM(i) = beta(i) / (omega(i) * epsilon_0 * epsilon_r);
        lambda_g(i) = 2 * pi / beta(i);
    end
end

% Plot Propagation Constant
figure;
plot(f / 1e9, beta, 'LineWidth', 2);
title('Propagation Constant vs Frequency');
xlabel('Frequency (GHz)');
ylabel('\beta (Propagation Constant, rad/m)');
grid on;

% Plot Group Velocity
figure;
plot(f / 1e9, v_g / 1e8, 'LineWidth', 2);
title('Group Velocity vs Frequency');
xlabel('Frequency (GHz)');
ylabel('v_g (Group Velocity, x10^8 m/s)');
grid on;

% Plot Phase Velocity
figure;
plot(f / 1e9, v_f / 1e8, 'LineWidth', 2);
title('Phase Velocity vs Frequency');
xlabel('Frequency (GHz)');
ylabel('v_f (Phase Velocity, x10^8 m/s)');
grid on;

% Plot Modal Impedance
figure;
plot(f / 1e9, Z_TM, 'LineWidth', 2);
title('Modal Impedance vs Frequency');
xlabel('Frequency (GHz)');
ylabel('Z_{TM} (Ohms)');
grid on;

% Plot Wavelength Inside Waveguide
figure;
plot(f / 1e9, lambda_g * 1e3, 'LineWidth', 2);
title('Wavelength Inside Waveguide vs Frequency');
xlabel('Frequency (GHz)');
ylabel('\lambda_g (Wavelength, mm)');
grid on;
