f_center = 1172e6; % Center frequency in Hz
freq = linspace(0.2*f_center, 1.8*f_center, 500); % Frequency range (Hz)
num_samples = 25; % Number of samples
variation_range = 0.06; % ±6%

% Initialize plots
figure;
hold on;

for sample = 1:num_samples
    % Random variations in Z_o1 and Z_o2
    Z_o1_varied = 42.43 * (1 + (rand()*2 - 1) * variation_range); % ±6% variation
    Z_o2_varied = 60 * (1 + (rand()*2 - 1) * variation_range);
    
    % Example varied S-Parameters
    S31_varied = -20 * log10(1 + freq / f_center * Z_o1_varied / 42.43);
    S32_varied = -20 * log10(1 + freq / f_center * 0.5 * Z_o2_varied / 60);
    
    % Plot the varied S31
    plot(freq/1e6, S31_varied, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5); % Gray lines for variations
end

% Add nominal curve for reference
S31_nominal = -20 * log10(1 + freq / f_center);
plot(freq/1e6, S31_nominal, 'r', 'LineWidth', 2); % Nominal case in red
xlabel('Frequency (MHz)');
ylabel('S31 (dB)');
title('S31 with Random Impedance Variations (±6%)');
grid on;
hold off;
