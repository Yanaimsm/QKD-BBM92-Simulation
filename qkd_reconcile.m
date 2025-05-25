% qkd_reconcile.m
% Step 3: Error correction and privacy amplification

% Load data generated in the previous step (raw keys and metadata)
load('qkd_coincidence.mat');

% --- QBER Estimation and Key Sifting ---

sample_size = floor(0.1 * valid_trials_count);      % Use 10% of raw key for QBER estimation
qber_sample_alice = alice_raw_key(1:sample_size);   % Slice of Alice's key
qber_sample_bob   = bob_raw_key(1:sample_size);     % Slice of Bob's key
errors = sum(qber_sample_alice ~= qber_sample_bob); % Count mismatches
qber = (errors / sample_size) * 100;                % Compute QBER as a percentage

% Display QBER to user
fprintf('Estimated QBER: %.2f%%\n', qber);

% Remove sampled portion from both keys (used for QBER testing)
alice_raw_key(1:sample_size) = [];
bob_raw_key(1:sample_size)   = [];
total_keys = length(alice_raw_key);  % Update usable key length

% Abort if QBER is too high (i.e., channel is insecure)
if qber > 11
    disp('Error: QBER exceeds threshold. Key distribution aborted.');
    disp(qber);
    return;
end

% --- Cascade Error Correction Setup ---

num_passes = 5;                          % Number of Cascade passes
k1 = 0.73 / (qber / 100);                % Initial block size (in bits), from Brassard-Salvail
block_sizes = ceil([k1, k1/2, k1/4, k1/8, k1/16]);  % Geometrically shrinking block sizes
alice_corrected_key = alice_raw_key;     % Initialize corrected key (copy)
bob_corrected_key   = bob_raw_key;       % Initialize corrected key (copy)
info_leakage = 0;                         % Tracks total number of leaked parity bits

% --- Cascade Error Correction ---

for pass = 1:num_passes
    % Shuffle both keys identically
    perm = randperm(total_keys);
    alice_corrected_key = alice_corrected_key(perm);
    bob_corrected_key   = bob_corrected_key(perm);

    block_size = block_sizes(pass);                % Block size for this pass
    num_blocks = floor(ceil(total_keys / block_size));

    for block_idx = 1:num_blocks
        start_idx = (block_idx - 1) * block_size + 1;
        end_idx   = min(block_idx * block_size, total_keys);
        if start_idx > total_keys || end_idx > total_keys
            continue;
        end

        % Slice current block
        a_block = alice_corrected_key(start_idx:end_idx);
        b_block = bob_corrected_key(start_idx:end_idx);

        % Compare parities of the two blocks
        if mod(sum(a_block), 2) ~= mod(sum(b_block), 2)
            info_leakage = info_leakage + 1;  % One parity leaked

            % Perform binary search to locate single error
            low = start_idx;
            high = end_idx;
            while low < high
                info_leakage = info_leakage + 1;  % One more parity leaked
                mid = floor((low + high)/2);
                if mod(sum(alice_corrected_key(low:mid)), 2) ~= mod(sum(bob_corrected_key(low:mid)), 2)
                    high = mid;
                else
                    low = mid + 1;
                end
            end

            % Flip the bit at the identified position
            bob_corrected_key(low) = ~bob_corrected_key(low);
        end
    end
end

% --- Privacy Amplification ---

n = length(alice_corrected_key);                % Final key length before PA
epsilon = 1e-9;                                  % Security parameter
s = floor(n - info_leakage - 2 * log2(1/epsilon));  % Length of final key after PA

% Sanity check: if final key is too short, abort
if s <= 0
    error('Privacy amplification failed: estimated final key length is non-positive.');
end

% Generate a random binary (s x n) matrix as universal hash
P = randi([0,1], s, n);

% Hash both Alice and Bob's corrected keys
KA = mod(P * alice_corrected_key, 2);
KB = mod(P * bob_corrected_key,   2);

% Verify that final keys match (they must!)
if any(KA ~= KB)
    error('Privacy amplification failed: mismatches remain.');
end

% Display final result

disp(['Final key length: ', num2str(s)]);
disp('Final key generation successful with privacy amplification.');
