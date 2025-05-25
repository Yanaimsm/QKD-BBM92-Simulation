% qkd_coincidence.m
% Step 2: Compute coincidences and extract raw key

load('qkd_params.mat');  % Load shared simulation parameters and variables

% --- Project amplitude for each basis combination (Alice/Bob: Z or X) ---

% Both use X basis
ampXX = PSY(1) * kron(Rot(alice_X_H)*H, Rot(bob_X_H)*H) + ...
        PSY(2) * kron(Rot(alice_X_H)*H, Rot(bob_X_H)*V) + ...
        PSY(3) * kron(Rot(alice_X_H)*V, Rot(bob_X_H)*H) + ...
        PSY(4) * kron(Rot(alice_X_H)*V, Rot(bob_X_H)*V);

% Alice uses X, Bob uses Z
ampXZ = PSY(1) * kron(Rot(alice_X_H)*H, Rot(bob_Z_H)*H) + ...
        PSY(2) * kron(Rot(alice_X_H)*H, Rot(bob_Z_H)*V) + ...
        PSY(3) * kron(Rot(alice_X_H)*V, Rot(bob_Z_H)*H) + ...
        PSY(4) * kron(Rot(alice_X_H)*V, Rot(bob_Z_H)*V);

% Alice uses Z, Bob uses X
ampZX = PSY(1) * kron(Rot(alice_Z_H)*H, Rot(bob_X_H)*H) + ...
        PSY(2) * kron(Rot(alice_Z_H)*H, Rot(bob_X_H)*V) + ...
        PSY(3) * kron(Rot(alice_Z_H)*V, Rot(bob_X_H)*H) + ...
        PSY(4) * kron(Rot(alice_Z_H)*V, Rot(bob_X_H)*V);

% Both use Z basis
ampZZ = PSY(1) * kron(Rot(alice_Z_H)*H, Rot(bob_Z_H)*H) + ...
        PSY(2) * kron(Rot(alice_Z_H)*H, Rot(bob_Z_H)*V) + ...
        PSY(3) * kron(Rot(alice_Z_H)*V, Rot(bob_Z_H)*H) + ...
        PSY(4) * kron(Rot(alice_Z_H)*V, Rot(bob_Z_H)*V);

% Convert amplitudes to detection probabilities
prob_coin_XX = abs(ampXX).^2;
prob_coin_XZ = abs(ampXZ).^2;
prob_coin_ZX = abs(ampZX).^2;
prob_coin_ZZ = abs(ampZZ).^2;

% Assign coincidence probabilities based on basis combinations
prob_coin = zeros(num_trials, 1);
prob_coin(alice_basis == 0 & bob_basis == 0) = prob_coin_ZZ(spd_rand(alice_basis == 0 & bob_basis == 0));
prob_coin(alice_basis == 0 & bob_basis == 1) = prob_coin_ZX(spd_rand(alice_basis == 0 & bob_basis == 1));
prob_coin(alice_basis == 1 & bob_basis == 0) = prob_coin_XZ(spd_rand(alice_basis == 1 & bob_basis == 0));
prob_coin(alice_basis == 1 & bob_basis == 1) = prob_coin_XX(spd_rand(alice_basis == 1 & bob_basis == 1));

% --- Determine output bits for Alice and Bob ---

% Default to 0
alice_outcome = zeros(size(spd_rand));
bob_outcome   = zeros(size(spd_rand));

% Based on SPD port number (1-4), determine outcomes
alice_outcome(spd_rand == 3 | spd_rand == 4) = 1;
bob_outcome(spd_rand == 2 | spd_rand == 4)   = 1;

% --- Coincidence detection (true, dark, and background) ---

% True coincidence: both detectors register detections
C = (random_coincidence <= prob_coin) & ...
    (random_detection_alice <= det_eff_alice) & ...
    (random_detection_bob   <= det_eff_bob);

% Dark noise: simultaneous false detections from detectors
dark_noise_alice = rand(num_trials, 1) < dark_rate_alice;
dark_noise_bob   = rand(num_trials, 1) < dark_rate_bob;
dark_coinc = dark_noise_alice & dark_noise_bob;

% Background noise: false counts due to environmental light
bg_coince = (rand_bg_noise_alice < bg_noise) & (rand_bg_noise_bob < bg_noise);

% Combine all coincidence sources
C = C | dark_coinc | bg_coince;

% --- Sifting: Keep only events where bases matched and a coincidence occurred ---

base_coin = (alice_basis == bob_basis);        % Basis agreement flag
valid_trials = (base_coin == 1) & (C == 1);     % Use only basis-matched, coincident events

% Extract raw keys from valid events
alice_raw_key = alice_outcome(valid_trials);
bob_raw_key   = bob_outcome(valid_trials);
valid_trials_count = length(alice_raw_key);

% Compute QBER on a subset (10%)
sample_size = floor(0.1 * valid_trials_count);
qber = (sum(alice_raw_key(1:sample_size) ~= bob_raw_key(1:sample_size)) / sample_size) * 100;

% Remove QBER test bits
alice_raw_key(1:sample_size) = [];
bob_raw_key(1:sample_size) = [];
total_keys = length(alice_raw_key);

% Compute recommended initial block size for Cascade error correction
block_size = ceil(0.73 / (qber / 100));  % Eq. from Brassard-Salvail

% Save results for next stage
save('qkd_coincidence.mat', ...
    'alice_raw_key', 'bob_raw_key', ...
    'valid_trials_count', 'qber', 'total_keys', 'block_size');
