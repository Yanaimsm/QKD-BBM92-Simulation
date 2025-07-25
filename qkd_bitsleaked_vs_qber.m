% QKD-BBM92 Simulation: Bits Leaked vs QBER (Cascade)
clear; clc; close all;

% PARAMETERS
N = 10000; % Number of raw key bits (keep as in article for comparison)
qber_list = 0.01:0.005:0.2; % QBER values (0.01 = 1%, 0.1 = 10%, ...)
num_trials = 100; % More trials per QBER value for better statistics
bits_leaked_results = zeros(length(qber_list), num_trials);

num_passes = 4;
% For storing error correction stats at QBER = 0.1
error_correction_qb010 = zeros(num_passes, num_trials); % [pass, trial]
total_errors_qb010 = zeros(1, num_trials); % For true error count
errors_remaining_per_qber = zeros(length(qber_list), num_trials);
for qi = 1:length(qber_list)
    qber = qber_list(qi);
    fprintf('Simulating QBER = %.3f\n', qber);
    for trial = 1:num_trials
        %% 1. Generate keys
        alice_key = randi([0 1], 1, N);
        bob_key   = alice_key;
        n_flip = round(qber * N);
        flip_idx = randperm(N, n_flip);
        bob_key(flip_idx) = ~bob_key(flip_idx);

        %% 2. Inline Cascade
        alice_raw_key = alice_key;
        bob_raw_key   = bob_key;

        % Simulate "QBER estimation sample"
        sample_size = floor(0.1 * N);
        alice_raw_key(1:sample_size) = [];
        bob_raw_key(1:sample_size)   = [];
        total_keys = length(alice_raw_key);

        % Cascade parameters
        k1 = 0.73 / (qber);
        block_sizes = ceil([k1, k1/2, k1*4, k1*8]);
        alice_corrected_key = alice_raw_key;
        bob_corrected_key   = bob_raw_key;
        info_leakage = 0;
        errors_corrected_pass = zeros(1, num_passes); % For error correction stat at QBER=0.1

        for pass = 1:num_passes
            perm = randperm(total_keys);
            alice_corrected_key = alice_corrected_key(perm);
            bob_corrected_key   = bob_corrected_key(perm);
            block_size = block_sizes(pass);
            num_blocks = floor(ceil(total_keys / block_size));
            for block_idx = 1:num_blocks
                start_idx = (block_idx - 1) * block_size + 1;
                end_idx   = min(block_idx * block_size, total_keys);
                if start_idx > total_keys || end_idx > total_keys
                    continue;
                end
                a_block = alice_corrected_key(start_idx:end_idx);
                b_block = bob_corrected_key(start_idx:end_idx);
                if mod(sum(a_block), 2) ~= mod(sum(b_block), 2)
                    info_leakage = info_leakage + 1;
                    errors_corrected_pass(pass) = errors_corrected_pass(pass) + 1;
                    low = start_idx;
                    high = end_idx;
                    while low < high
                        info_leakage = info_leakage + 1;
                        mid = floor((low + high)/2);
                        if mod(sum(alice_corrected_key(low:mid)), 2) ~= mod(sum(bob_corrected_key(low:mid)), 2)
                            high = mid;
                        else
                            low = mid + 1;
                        end
                    end
                    bob_corrected_key(low) = ~bob_corrected_key(low);
                end
            end
        end

        if abs(qber - 0.1) < 1e-6
            error_correction_qb010(:, trial) = errors_corrected_pass(:);
            total_errors_qb010(trial) = n_flip;
        end

        % Store the number of bits leaked in this trial
        bits_leaked_results(qi, trial) = info_leakage;

        % Store errors remaining in the key after all passes
        errors_remaining_per_qber(qi, trial) = sum(alice_corrected_key ~= bob_corrected_key);
    end
end

% --- Print error correction stats for QBER = 0.1 ---
if any(abs(qber_list - 0.1) < 1e-6)
    mean_corr = mean(error_correction_qb010, 2);
    mean_total_errors = mean(total_errors_qb010);
    fprintf('\n--- Error Correction for QBER = 0.1 ---\n');
    for pass = 1:num_passes
        perc = 100 * mean_corr(pass) / mean_total_errors;
        fprintf('Pass %d: %.2f%% errors corrected\n', pass, perc);
    end
    fprintf('Total errors per trial: %.2f\n', mean_total_errors);
end

% --- Print errors remaining per QBER:
fprintf('\n--- Average errors remaining after Cascade per QBER ---\n');
for qi = 1:length(qber_list)
    avg_err_remain = mean(errors_remaining_per_qber(qi, :));
    fprintf('QBER = %.3f: Average errors remaining = %.2f\n', qber_list(qi), avg_err_remain);
end

% --- After your main simulation loop ---
bits_leaked_mean = mean(bits_leaked_results, 2);
bits_leaked_std = std(bits_leaked_results, 0, 2);

% Calculate the upper bound D (from the formula) for each QBER:
D_bound = zeros(size(qber_list));
N_sifted = N - floor(0.1 * N); % After QBER estimation

for qi = 1:length(qber_list)
    qber = qber_list(qi);
    k1 = 0.73 / (qber);
    block_sizes = [k1, k1*2, k1*4, k1*8];
    D = 0;
    for pass = 1:4
        k_i = block_sizes(pass);
        num_blocks = N_sifted / k_i;
        N_errors = qber * N_sifted / 4;
        D = D + num_blocks + N_errors * log2(k_i);
    end
    D_bound(qi) = D;
end

qber_list_for_plot = qber_list * 100; % percent for plotting

% --- Plot ---
figure;
errorbar(qber_list_for_plot, bits_leaked_mean, bits_leaked_std, 'r', 'LineWidth', 1.5); hold on;
plot(qber_list_for_plot, D_bound, 'k--', 'LineWidth', 1.5);
xlabel('QBER (%)');
ylabel('Bits leaked');
title('Bits leakage in Cascade error correction algorithm for QKD');
grid on;
legend({'Cascade (Sim)', 'Cascade theoretical upper bound (Eq.)'}, 'Location', 'northwest');

% Show values at QBER = 1%, 5%, 11%, 17% as (x, y), all labels in black and offset
qber_marks = [0.01, 0.05, 0.11, 0.17] * 100; % percent
y_offset = [-300, -200, 200, -200]; % For better spacing

for i = 1:length(qber_marks)
    [~, idx] = min(abs(qber_list_for_plot - qber_marks(i)));
    sim_val = bits_leaked_mean(idx);
    bound_val = D_bound(idx);
    plot(qber_list_for_plot(idx), sim_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'HandleVisibility','off');
    text(qber_list_for_plot(idx)+0.5, sim_val + y_offset(i), ...
        sprintf('(%.0f, %.0f)', qber_list_for_plot(idx), sim_val), ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
    plot(qber_list_for_plot(idx), bound_val, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'HandleVisibility','off');
    text(qber_list_for_plot(idx)+0.5, bound_val - y_offset(i), ...
        sprintf('(%.0f, %.0f)', qber_list_for_plot(idx), bound_val), ...
        'Color', 'k', 'FontSize', 10, 'FontWeight', 'bold');
end
