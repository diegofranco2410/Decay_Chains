%% MOLYBDENUM-99 DECAY CHAIN SIMULATION (Medical Isotope Generator)
%  This script simulates the decay kinetics of the Mo-99 -> Tc-99m generator.
%  It highlights the "Transient Equilibrium" phenomenon and the branching
%  ratio between the metastable state (Tc-99m) and the ground state (Tc-99).
%
%  Chain:
%  1. Mo-99  (Parent, T1/2 = 65.9 h)  -> Beta decay
%  2. Tc-99m (Meta,   T1/2 = 6.01 h)  -> Isomeric Transition (Gamma)
%  3. Tc-99  (Ground, T1/2 = 2.1e5 y) -> Beta decay (Quasi-stable)
%  4. Ru-99  (Stable)

clear; clc; close all;

% 1. Physics Setup
names = {'Mo-99', 'Tc-99m', 'Tc-99 (Ground)', 'Ru-99 (Stable)'};

% Time units conversion (to Hours)
h = 1;
y = 365.25 * 24; 

% Half-lives (in Hours)
% Note: Tc-99 ground state is essentially stable for this timeframe
half_lives = [65.94*h, 6.01*h, 2.111e5*y]; 

% Initial Activity Calculation
% In medical physics, we usually start with Activity (Ci/Bq), not N.
% Let's assume a fresh 15 Curie (15 Ci) generator.
A0_Ci = 15; 
lambda_Mo = log(2) / half_lives(1);
% Convert Activity to Number of atoms: N = A / lambda
% (1 Ci = 3.7e10 disintegrations/s). We need lambda in seconds for this conversion.
N0 = (A0_Ci * 3.7e10) / (lambda_Mo / 3600); 

fprintf('Simulating a %.1f Ci Mo-99 Generator...\n', A0_Ci);
fprintf('Initial Mo-99 Atoms: %.2e\n', N0);

% --- TOPOLOGY (Branching) ---
% Matrix size: 3 Parents -> 4 Total Species
B = zeros(3, 4);

% Branch 1: Mo-99 (Idx 1)
% ~87.6% goes to metastable Tc-99m (Idx 2)
% ~12.4% goes directly to ground Tc-99 (Idx 3)
B(1, 2) = 0.876; 
B(1, 3) = 0.124;

% Branch 2: Tc-99m (Idx 2)
% 100% decays to ground state Tc-99 (Idx 3)
B(2, 3) = 1.0;

% Branch 3: Tc-99 Ground (Idx 3)
% 100% decays to Ru-99 (Idx 4)
B(3, 4) = 1.0;

% 2. Time Configuration
% We simulate 1 week (168 hours) to see the peak and equilibrium
T = 0:0.01:168; % Linear scale (Hours)

% 3. Run Models
% Analytical Solution
[N, A, ~] = DecayMatrix(N0, T, half_lives, ...
    'Branching', B, ...
    'Nuclei', names, ...
    'TimeUnits', 'h', ... % Important for correct Activity calculation
    'Plotting', false);

% Monte Carlo (Scaled down for visibility of fluctuations)
% Since N0 is ~10^16, MC is too slow/heavy for exact N0. 
% We scale down N0 for MC to show the statistical nature, 
% or we skip MC if we want to show exact medical curves.
% Let's do a small MC run just to show the logic works.
N0_mc = 5000; 
[N_mc, A_mc, ~] = DecayMonteCarlo(N0_mc, T, half_lives, ...
    'Branching', B, ...
    'Nuclei', names, ...
    'TimeUnits', 'h', ...
    'Plotting', false);

% 4. Visualization 
figure('Color', 'w', 'Position', [100 100 1000 800]);

% --- Subplot 1: ACTIVITY (The most important for MedPhys) ---
subplot(2,1,1);
hold on; grid on; box on;

% Plot Analytical Activities
% Mo-99 (Parent)
plot(T, A(1,:), 'b-', 'LineWidth', 2); 
% Tc-99m (Daughter - The useful one)
plot(T, A(2,:), 'r-', 'LineWidth', 2); 
% Tc-99 (Ground - Useless waste)
plot(T, A(3,:), 'g--', 'LineWidth', 1.5);

xlabel('Time (Hours)', 'FontSize', 12);
ylabel('Activity (Ci)', 'FontSize', 12);
title('Mo-99 / Tc-99m Generator Kinetics (Transient Equilibrium)', 'FontSize', 14);
legend({'Mo-99 (Parent)', 'Tc-99m (Metastable)', 'Tc-99 (Ground)'}, 'Location', 'best');

% Highlight the peak
[max_act, idx_max] = max(A(2,:));
t_max = T(idx_max);
xline(t_max, '--k', sprintf('Max Tc-99m yield\n(t = %.2f)', t_max), ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% Subplot 2: Monte Carlo Dynamics (Scaled) ---
subplot(2,1,2);
hold on; grid on; box on;
stairs(T, A_mc(1,:), 'b', 'LineWidth', 1);
stairs(T, A_mc(2,:), 'r', 'LineWidth', 1);
stairs(T, A_mc(3,:), 'g', 'LineWidth', 1);

xlabel('Time (Hours)', 'FontSize', 12);
ylabel('Activity [Scaled MC]', 'FontSize', 12);
title(['Stochastic Simulation (N_0 = ' num2str(N0_mc) ')'], 'FontSize', 14);
legend(names{1:3});

% 5. Console Analysis
fprintf('\n--- Generator Analysis ---\n');
fprintf('Theoretical max Tc-99m activity occurs at: %.2f hours\n', t_max);
fprintf('At peak, Tc-99m activity is %.2f%% of Mo-99 initial activity.\n', (max_act/A0_Ci)*100);