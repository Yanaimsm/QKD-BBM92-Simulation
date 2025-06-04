% Step 1: Parameters initialization

T_sim = 10;                      % Total simulation time in seconds
source_rate = 12e4;              % Photon pair generation rate (pairs/sec) from SPDC source
num_trials = round(source_rate * T_sim);  % Total number of photon pair trials to simulate

% Detection efficiencies (including optical losses)
det_eff_alice = 0.453;           % Overall detection probability at Alice
det_eff_bob = 0.167;             % Overall detection probability at Bob

% Detector dark count and environmental background noise
dark_rate_alice = 0.0463;        % Probability of false click due to dark noise at Alice
dark_rate_bob   = 0.0283;        % Probability of false click due to dark noise at Bob
bg_noise = 0.01;                 % Background light contribution probability (applied later)

% Define standard polarization states for qubit encoding
H = [1; 0];                      % Horizontal polarization
V = [0; 1];                      % Vertical polarization

% Define the entangled Bell state |Ψ⟩ = (|HH⟩ + |VV⟩)/√2
PSY = [1; 0; 0; 1] / sqrt(2);    % Tensor product representation: kron(H,H) + kron(V,V)

% Define a general 2D rotation matrix for rotating polarization bases
Rot = @(theta) [cosd(theta), -sind(theta); sind(theta), cosd(theta)];

% Define angles for Z and X basis measurements at Alice
alice_Z_H = 0;                   % Horizontal (0°) in Z basis
alice_Z_V = 90;                  % Vertical (90°) in Z basis
alice_X_H = 45;                  % +45° for X basis (|+⟩)
alice_X_V = 135;                 % −45° for X basis (|−⟩)

% Define angles for Z and X basis measurements at Bob (same as Alice)
bob_Z_H = 0;
bob_Z_V = 90;
bob_X_H = 45;
bob_X_V = 135;

% Random basis selection: 0 for Z basis, 1 for X basis
alice_basis = randi([0, 1], num_trials, 1);  % Alice randomly chooses Z or X
bob_basis   = randi([0, 1], num_trials, 1);  % Bob randomly chooses Z or X

% Simulated photodetector response: random assignment to SPD channel 1 to 4
spd_rand = randi([1, 4], num_trials, 1);     % Used to resolve click outcomes

% Generate random values to simulate noise, detection probability, and coincidence logic
random_coincidence      = rand(num_trials, 1);  % Probability of SPDC-based coincidence
random_detection_alice  = rand(num_trials, 1);  % Whether Alice's photon is detected
random_detection_bob    = rand(num_trials, 1);  % Whether Bob's photon is detected
rand_bg_noise_bob       = rand(num_trials, 1);  % Background noise events at Bob
rand_bg_noise_alice     = rand(num_trials, 1);  % Background noise events at Alice

% Save all relevant parameters into a .mat file for later use in other scripts
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
