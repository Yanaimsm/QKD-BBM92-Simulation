clc;
fprintf('QKD Simulation - Key Length & QBER vs Dark Count Rates (Logarithmic Scale)\n');
fprintf('---------------------------------------------------------------------------\n');

% --- Simulation Parameters ---
num_sweep_points = 20; % Run over 20 values
num_simulations_per_point = 5; % Simulate each value 5 times for averaging

% Define sweep range for dark rates using a logarithmic scale
% Going from 10^-5 (0.00001) to 10^-1 (0.1) with num_sweep_points
dark_rate_vals = logspace(-5, -1, num_sweep_points);

% Initialize storage for average results
avg_key_lengths = zeros(size(dark_rate_vals));
avg_qber_vals = zeros(size(dark_rate_vals));

% --- Outer Simulation Loop (over dark_rate_vals) ---
for i = 1:length(dark_rate_vals)
    current_dark_rate = dark_rate_vals(i);
    
    % Initialize temporary storage for results of inner simulations
    temp_key_lengths = zeros(1, num_simulations_per_point);
    temp_qber_vals = zeros(1, num_simulations_per_point);
    
    fprintf('\nSimulating Dark Rate: %.5f (Iteration %d/%d)\n', current_dark_rate, i, length(dark_rate_vals));
    
    % --- Inner Simulation Loop (for averaging) ---
    for j = 1:num_simulations_per_point
        fprintf('  Sub-simulation %d/%d...\n', j, num_simulations_per_point);
        
        % Load base parameters (this resets parameters for each sub-simulation)
        % Note: qkd_params.m is provided in the prompt. We assume it properly
        % generates the random variables (alice_basis, random_detection_alice, etc.)
        % each time it's run, which is crucial for statistical averaging.
        run('qkd_params.m');
        
        % Set current dark rates (same for Alice and Bob)
        dark_rate_alice = current_dark_rate;
        dark_rate_bob = current_dark_rate;
        
        % Save updated parameters to qkd_params.mat
        % This part is crucial as qkd_coincidence.m will load from this file.
        save('qkd_params.mat', ...
            'T_sim', 'source_rate', 'num_trials', ...
            'det_eff_alice', 'det_eff_bob', ...
            'dark_rate_alice', 'dark_rate_bob', ...
            'bg_noise', ...
            'H', 'V', 'PSY', ...
            'Rot', ...
            'alice_Z_H', 'alice_Z_V', 'alice_X_H', 'alice_X_V', ...
            'bob_Z_H', 'bob_Z_V', 'bob_X_H', 'bob_X_V', ...
            'alice_basis', 'bob_basis', 'spd_rand', ...
            'random_coincidence', 'random_detection_alice', 'random_detection_bob', ...
            'rand_bg_noise_bob', 'rand_bg_noise_alice');
        
        % Run coincidence detection
        run('qkd_coincidence.m');
        
        try
            % Run reconciliation (qkd_reconcile.m must calculate 'qber' and 's')
            run('qkd_reconcile.m');
            
            % Check if 's' is defined and has a valid value
            if exist('s', 'var') && s > 0
                temp_key_lengths(j) = s;
            else
                % If s is 0 or not defined, set key length to 0 for this sub-simulation
                temp_key_lengths(j) = 0;
                warning('qkd_reconcile.m did not define a positive ''s'' value or s was 0 for this sub-simulation.');
            end
            
            % Ensure 'qber' is available from qkd_reconcile.m
            if exist('qber', 'var')
                temp_qber_vals(j) = qber;
            else
                warning('qkd_reconcile.m did not define ''qber'' for this sub-simulation.');
                temp_qber_vals(j) = NaN;
            end
            
        catch ME
            fprintf('  Sub-simulation failed: %s\n', ME.message);
            temp_key_lengths(j) = 0; % Set to 0 if an error occurs
            temp_qber_vals(j) = NaN; % Set to NaN if an error occurs
        end
    end % End of inner simulation loop
    
    % --- Calculate averages for the current dark_rate_vals(i) ---
    % Filter out zeros/NaNs before averaging for key length (only average successful key generations)
    valid_keys = temp_key_lengths(temp_key_lengths > 0);
    if ~isempty(valid_keys)
        avg_key_lengths(i) = mean(valid_keys);
    else
        avg_key_lengths(i) = 0; % No successful key generations for this dark rate
    end

    % Filter out NaNs for QBER (only average successful QBER calculations)
    valid_qbers = temp_qber_vals(~isnan(temp_qber_vals));
    if ~isempty(valid_qbers)
        avg_qber_vals(i) = mean(valid_qbers);
    else
        avg_qber_vals(i) = NaN; % No successful QBER calculations for this dark rate
    end
    
    fprintf('  Average Key length: %d | Average QBER: %.2f%%\n', ...
            round(avg_key_lengths(i)), avg_qber_vals(i)); % Round key length for display
end % End of outer simulation loop

% Create Combined Plot

figure;

% Filter out points where avg_key_lengths is 0 for key length plot
valid_indices_key = avg_key_lengths > 0;
% Filter out NaN QBER values for QBER plot
valid_indices_qber = ~isnan(avg_qber_vals);

% Plot Key Length on the left y-axis
yyaxis left;
semilogx(dark_rate_vals(valid_indices_key), avg_key_lengths(valid_indices_key), 'b-o', 'LineWidth', 2, 'DisplayName', 'Average Final Key Length');
ylabel('Average Final Key Length (bits)');
hold on; % Keep the plot active to add the second curve

% Plot QBER on the right y-axis
yyaxis right;
semilogx(dark_rate_vals(valid_indices_qber), avg_qber_vals(valid_indices_qber), 'r-x', 'LineWidth', 2, 'DisplayName', 'Average QBER');
ylabel('Average QBER (%)');

% --- Customize X-axis for Percentage Display ---
ax = gca; % Get current axes handle (applies to both y-axes)
xticks_loc = get(ax, 'XTick');
% Convert tick values to percentage strings
xticks_labels = arrayfun(@(x) sprintf('%.3g%%', x * 100), xticks_loc, 'UniformOutput', false);
set(ax, 'XTickLabel', xticks_labels);

% Set common X-axis label and title
xlabel('Dark Count Rate (Alice and Bob)');
title('Average Key Length and QBER vs. Dark Count Rate');
grid on;
legend('Location', 'best'); % Display legend for both curves

% Adjust figure position
set(gcf, 'Position', [100 100 800 600]);