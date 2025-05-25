% qkd_simulate.m – Run simulation sweep based on user choice

clc;
fprintf('QKD Parameter Sweep Simulation\n');
fprintf('------------------------------\n');
fprintf('Select parameter to sweep:\n');
fprintf('1. Background Noise (bg_noise)\n');
fprintf('2. Alice Detection Efficiency (det_eff_alice)\n');
fprintf('3. Alice Dark Count Rate (dark_rate_alice)\n');
choice = input('Enter your choice (1/2/3): ');

% Define sweep variable
switch choice
    case 1
        param_name = 'bg_noise';
        sweep_vals = linspace(0, 0.1, 20);
    case 2
        param_name = 'det_eff_alice';
        sweep_vals = linspace(0.1, 0.9, 20);
    case 3
        param_name = 'dark_rate_alice';
        sweep_vals = linspace(0, 0.1, 20);
    otherwise
        error('Invalid selection. Please enter 1, 2, or 3.');
end

% Preallocate QBER storage
qber_vals = NaN(size(sweep_vals));

% Sweep simulation loop
for i = 1:length(sweep_vals)
    run('qkd_params.m');        % Load base parameters
    eval([param_name ' = sweep_vals(i);']);  % Override selected parameter
    save('qkd_params.mat');     % Save updated parameter
    
    run('qkd_coincidence.m');   % Generate coincidences and raw keys
    try
        run('qkd_reconcile.m');     % Run error correction and QBER calc
        qber_vals(i) = qber;        % Store result
    catch
        fprintf('Failed at %s = %.4f\n', param_name, sweep_vals(i));
    end
    fprintf('current %s = %.4f\n', param_name, sweep_vals(i));
end

% Plot QBER vs Selected Parameter
figure;
plot(sweep_vals, qber_vals, '-o', 'LineWidth', 2);
xlabel(strrep(param_name, '_', '\_'));
ylabel('QBER (%)');
title(sprintf('QBER vs. %s', strrep(param_name, '_', '\_')));
grid on;
