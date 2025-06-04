% qkd_display.m – Master QKD script to run simulation and show results
clc; close all;
disp('Quantum Key Distribution Simulation');
disp('===================================');

% Step 1: Initialize parameters
run('qkd_params.m');

% Step 2: Simulate photon pair detection and generate raw keys
run('qkd_coincidence.m');

% Step 3: Run error correction and privacy amplification
run('qkd_reconcile.m');

% --- Step 4: Display results in a tabbed figure ---

% Create a figure
hFig = figure('Name', 'QKD Simulation Results', 'NumberTitle', 'off', 'Position', [100, 100, 900, 700]);

% Create a tab group
hTabGroup = uitabgroup(hFig);

%% Tab 1: Simulation Summary
hTab1 = uitab(hTabGroup, 'Title', 'Simulation Summary');
hPanel1 = uipanel(hTab1, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Display summary information using uicontrols (text)
uicontrol(hPanel1, 'Style', 'text', 'String', '=== Simulation Summary ===', ...
    'Units', 'normalized', 'Position', [0.05 0.9 0.9 0.05], 'HorizontalAlignment', 'left', 'FontSize', 12, 'FontWeight', 'bold');

summary_str = sprintf([...
    'Simulation time:        %.2f seconds\n', ...
    'Total trials:           %d\n', ...
    'Valid coincidence trials: %d\n', ...
    'Raw key length:         %d bits\n', ...
    'Estimated QBER:         %.2f %%\n', ...
    'Final key length:       %d bits\n', ...
    'Final key accuracy:     %.2f %%\n', ...
    'Information leakage:    %d bits\n'], ...
    T_sim, num_trials, valid_trials_count, total_keys, qber, s, final_key_accuracy, info_leakage);

uicontrol(hPanel1, 'Style', 'text', 'String', summary_str, ...
    'Units', 'normalized', 'Position', [0.05 0.3 0.9 0.55], 'HorizontalAlignment', 'left', 'FontSize', 10, 'FontName', 'FixedWidth');

%% Tab 2: QBER Histogram
hTab2 = uitab(hTabGroup, 'Title', 'QBER Histogram');
hAxes2 = axes('Parent', hTab2, 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

% Plot histogram of mismatches on this tab's axes
bar(hAxes2, index_col, mismatch_col, 'FaceColor', [0.7 0 0]);
xlabel(hAxes2, 'Key Bit Index');
ylabel(hAxes2, 'Mismatch (1=True)');
title(hAxes2, 'Bit Mismatches After Privacy Amplification');
grid(hAxes2, 'on');

%% Tab 3: Final Key Data
hTab3 = uitab(hTabGroup, 'Title', 'Final Key Data');

% Prepare data for the uitable
% Show Alice's corrected key, Bob's corrected key, and their mismatches
column_headers = {'Index', 'Alice''s Final Key', 'Bob''s Final Key', 'Mismatch'};
% Demonstration, assuming alice_final_key_col and bob_final_key_col exist.
final_key_cell = cell(min_len, 4);
for i = 1:min_len
    final_key_cell{i, 1} = index_col(i);
    final_key_cell{i, 2} = alice_final_key_col(i);
    final_key_cell{i, 3} = bob_final_key_col(i);
    final_key_cell{i, 4} = mismatch_col(i);
end

uitable(hTab3, 'Data', final_key_cell, ...
        'ColumnName', column_headers, ...
        'RowName', {}, ...
        'Units', 'Normalized', ...
        'Position', [0, 0, 1, 1]);