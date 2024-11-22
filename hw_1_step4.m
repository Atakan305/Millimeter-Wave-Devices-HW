% Constants
c = 3e8; % Speed of light (m/s)
a = 0.0740; % Waveguide width (m)
b = 0.0296; % Waveguide height (m)
epsilon_r = 3.7; % Relative permittivity
epsilon_0 = 8.854e-12; % Vacuum permittivity (F/m)
mu_0 = 4 * pi * 1e-7; % Vacuum permeability (H/m)
f_c = (c / 2) * sqrt((3 / a)^2 + (4 / b)^2) / sqrt(epsilon_r); % Cutoff frequency
f = 1.5 * f_c; % Operating frequency (1.5 times cutoff frequency)
omega = 2 * pi * f;
k_c = sqrt((3 * pi / a)^2 + (4 * pi / b)^2); % Cutoff wavenumber
k = omega * sqrt(mu_0 * epsilon_0 * epsilon_r);
beta = sqrt(k^2 - k_c^2); % Propagation constant

% Spatial grid
x = linspace(0, a, 100);
y = linspace(0, b, 100);
[X, Y] = meshgrid(x, y);
z = 0; % Assume z = 0 for 2D field plotting

% Longitudinal Electric Field (E_z)
E0 = 1; % Amplitude
Ez = E0 * sin(3 * pi * X / a) .* sin(4 * pi * Y / b) .* exp(-1j * beta * z);

% Transverse Electric Fields (Ex, Ey)
Ex = -1j * beta * (3 * pi / a) * E0 .* cos(3 * pi * X / a) .* sin(4 * pi * Y / b) ./ k_c^2 .* exp(-1j * beta * z);
Ey = -1j * beta * (4 * pi / b) * E0 .* sin(3 * pi * X / a) .* cos(4 * pi * Y / b) ./ k_c^2 .* exp(-1j * beta * z);

% Magnetic Fields (Hx, Hy)
Hx = (1j * omega * epsilon_0 * epsilon_r * (4 * pi / b)) * E0 .* sin(3 * pi * X / a) .* cos(4 * pi * Y / b) ./ k_c^2 .* exp(-1j * beta * z);
Hy = -(1j * omega * epsilon_0 * epsilon_r * (3 * pi / a)) * E0 .* cos(3 * pi * X / a) .* sin(4 * pi * Y / b) ./ k_c^2 .* exp(-1j * beta * z);

% Active Power (Poynting Vector, S = E x H*)
Sx = real(Ez .* conj(Hy));
Sy = real(-Ez .* conj(Hx));

% Debugging: Print max values for verification
disp(['Max |Ex|: ', num2str(max(max(abs(real(Ex)))))]);
disp(['Max |Ey|: ', num2str(max(max(abs(real(Ey)))))]);
disp(['Max |Hx|: ', num2str(max(max(abs(real(Hx)))))]);
disp(['Max |Hy|: ', num2str(max(max(abs(real(Hy)))))]);
disp(['Max |Sx|: ', num2str(max(max(abs(Sx))))]);
disp(['Max |Sy|: ', num2str(max(max(abs(Sy))))]);

% Apply scaling factor for visualization
scaling_factor = 1e6;
Ex_scaled = real(Ex) * scaling_factor;
Ey_scaled = real(Ey) * scaling_factor;
Hx_scaled = real(Hx) * scaling_factor;
Hy_scaled = real(Hy) * scaling_factor;
Sx_scaled = Sx * scaling_factor;
Sy_scaled = Sy * scaling_factor;

% Reduce quiver density
skip = 5;

% Plot Longitudinal Electric Field (E_z)
figure;
contourf(X, Y, real(Ez), 50, 'LineStyle', 'none');
colorbar;
title('Longitudinal Electric Field (E_z) for TM_{3,4} Mode');
xlabel('x (m)');
ylabel('y (m)');

% Plot Transverse Electric Fields (Ex, Ey)
figure;
quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), Ex_scaled(1:skip:end, 1:skip:end), Ey_scaled(1:skip:end, 1:skip:end), 'b');
title('Transverse Electric Field (Ex, Ey) for TM_{3,4} Mode');
xlabel('x (m)');
ylabel('y (m)');

% Plot Transverse Magnetic Fields (Hx, Hy)
figure;
quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), Hx_scaled(1:skip:end, 1:skip:end), Hy_scaled(1:skip:end, 1:skip:end), 'r');
title('Transverse Magnetic Field (Hx, Hy) for TM_{3,4} Mode');
xlabel('x (m)');
ylabel('y (m)');

% Plot Active Power (Poynting Vector)
figure;
quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), Sx_scaled(1:skip:end, 1:skip:end), Sy_scaled(1:skip:end, 1:skip:end), 'k');
title('Active Power Flow (Poynting Vector) for TM_{3,4} Mode');
xlabel('x (m)');
ylabel('y (m)');
