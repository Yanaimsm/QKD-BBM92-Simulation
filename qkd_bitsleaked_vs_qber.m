% QKD-BBM92 Simulation: Bits Leaked vs QBER (Cascade)
clear; clc; close all;

% PARAMETERS
N = 10000; % Number of raw key bits (keep as in article for comparison)
qber_list = 0.01:0.005:0.17; % QBER values: more points, up to 0.22 (where Cascade typically fails)
num_trials = 50; % More trials per QBER value for better statistics
bits_leaked_results = zeros(length(qber_list), num_trials);


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

        %% 2. Inline Cascade (from your script)
        % Use all of Alice and Bob's keys directly
        alice_raw_key = alice_key;
        bob_raw_key   = bob_key;
        valid_trials_count = N;

        % Simulate "QBER estimation sample" (like your script)
        sample_size = floor(0.1 * valid_trials_count);
        alice_raw_key(1:sample_size) = [];
        bob_raw_key(1:sample_size)   = [];
        total_keys = length(alice_raw_key);

        % Cascade parameters
        num_passes = 4;
        %qber_for_block = qber * 100; % As percent
        k1 = 0.73 / (qber);
        block_sizes = ceil([k1, k1*2, k1*4, k1*8]);
        alice_corrected_key = alice_raw_key;
        bob_corrected_key   = bob_raw_key;
        info_leakage = 0;

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

        % Store the number of bits leaked in this trial
        bits_leaked_results(qi, trial) = info_leakage;
    end
end

% ANALYSIS
bits_leaked_mean = mean(bits_leaked_results, 2);
bits_leaked_std = std(bits_leaked_results, 0, 2);

% PLOT
figure;
errorbar(qber_list, bits_leaked_mean, bits_leaked_std, 'r', 'LineWidth', 1.5);
xlabel('QBER'); % QBER in fraction (0.xx)

ylabel('Bits leaked');
title('Error Reconciliation in QKD (Cascade)');
grid on;
legend({'Cascade'});
