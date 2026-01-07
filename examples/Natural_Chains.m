%% NATURAL DECAY CHAINS SIMULATION (Comparison: Analytic vs Monte Carlo)
%  This script simulates the 4 natural radioactive series:
%  1. Thorium Series (4n)
%  2. Neptunium Series (4n + 1)
%  3. Uranium Series (4n + 2)
%  4. Actinium Series (4n + 3)
%
%  It handles the complex branching ratios (e.g., Bi-212) and spans time
%  scales from microseconds to billions of years using logarithmic steps.

clear; clc; close all;

% 1. Select Chain
fprintf('Select a Natural Decay Chain to simulate:\n');
fprintf('1. Thorium Series (Th-232)\n');
fprintf('2. Neptunium Series (Np-237)\n');
fprintf('3. Uranium Series (U-238)\n');
fprintf('4. Actinium Series (U-235)\n');
choice = input('Enter number [1-4]: ');

[names, half_lives_s, B, description] = get_chain_data(choice);

% 2. Simulation Setup
% Physics Parameters
N0 = 1000; % Small number to visualize Monte Carlo fluctuations (Stochastic noise)
% If you want smooth curves, increase to 1e5 or 1e6.

% Time Configuration (Logarithmic Scale)
% We span from 1e-6 years to 1e11 years to see the whole history
years_to_sec = 3.154e7;
T_years = logspace(-8, 11, 1000); 
T = T_years * years_to_sec; % Convert to seconds for calculation

% 3. Run Models
fprintf('Simulating %s...\n', description);

% --- FIX NUMERICAL STIFFNESS ---
% Extremely short half-lives (ms) break numeric precision
% in geological time scales. Las limitamos a un minimo razonable (ej. 1 minuto).
min_stable_limit = 3600*5; % 5 hours in seconds
half_lives_s(half_lives_s < min_stable_limit) = min_stable_limit;
% -------------------------------

% A) Analytical Solution (Matrix Exponential)
fprintf('-> Calculating Analytical Solution...\n');
[N_ana, ~, ~] = DecayMatrix(N0, T, half_lives_s, ...
    'Branching', B, ...
    'Nuclei', names, ...
    'Plotting', false);

% B) Stochastic Solution (Monte Carlo)
fprintf('-> Calculating Monte Carlo Simulation...\n');
[N_mc, ~, ~] = DecayMonteCarlo(N0, T, half_lives_s, ...
    'Branching', B, ...
    'Nuclei', names, ...
    'Plotting', false);

% 4. Visualization
fig = figure('Color', 'w', 'Position', [50 50 1200 900]);
sgtitle(['Natural Decay Series: ' description], 'FontSize', 18, 'FontWeight', 'bold');

% --- Subplot 1: Population Dynamics ---
subplot(3,1,[1 2]); % Use 2/3 of the screen
hold on; box on; grid on;

% Color map
cmap = jet(length(names)); 

% Plot lines
pl = [];
for i = 1:length(names)
    % Monte Carlo (Thin lines / Stairs)
    semilogx(T_years, N_mc(i,:), '-', 'Color', [cmap(i,:) 0.4], 'LineWidth', 1);
    
    % Analytical (Thick lines)
    p = semilogx(T_years, N_ana(i,:), '-', 'Color', cmap(i,:), 'LineWidth', 2.5);
    pl = [pl, p];
end

set(gca, 'XScale', 'log');
ylabel('Population N(t)', 'FontSize', 12);
title(['Evolution from Microseconds to Billions of Years (N_0 = ' num2str(N0) ')'], 'FontSize', 14);
xlim([T_years(1) T_years(end)]);
ylim([0 N0*1.1]);

% Legend (Filter slightly if too many elements, or show all)
legend(pl, names, 'Location', 'eastoutside', 'FontSize', 8);

% --- Subplot 2: Residuals (Stochastic Noise) ---
subplot(3,1,3);
hold on; box on; grid on;

% Calculate Residuals (MC - Theory)
Residuals = N_mc - N_ana;

for i = 1:length(names)
    semilogx(T_years, Residuals(i,:), '.', 'Color', cmap(i,:), 'MarkerSize', 5);
end
yline(0, 'k-', 'LineWidth', 1.5);

set(gca, 'XScale', 'log');
xlabel('Time (Years) [Log Scale]', 'FontSize', 14);
ylabel('Residuals (\Delta N)', 'FontSize', 12);
title('Monte Carlo Fluctuations (Quantum Noise)', 'FontSize', 12);
xlim([T_years(1) T_years(end)]);

fprintf('Done.\n');

%
%  HELPER FUNCTION: DATABASE OF NATURAL CHAINS
function [names, half_lives, B, desc] = get_chain_data(selection)

    % Time Conversions
    y = 3.154e7;  % Year
    d = 86400;    % Day
    h = 3600;     % Hour
    m = 60;       % Minute
    s = 1;        % Second
    ms = 1e-3;    % Millisecond
    us = 1e-6;    % Microsecond
    ns = 1e-9;    % Nanosecond

    switch selection
        case 1 % THORIUM SERIES (4n)
            desc = 'Thorium Series (Th-232)';
            % Simplified dominant path + Bi-212 Branching
            names = {'Th-232', 'Ra-228', 'Ac-228', 'Th-228', 'Ra-224', ...
                     'Rn-220', 'Po-216', 'Pb-212', 'Bi-212', 'Po-212', 'Tl-208', 'Pb-208 (Stable)'};
            
            half_lives = [1.405e10*y, 5.75*y, 6.25*h, 1.9116*y, 3.6319*d, ...
                          55.6*s, 0.145*s, 10.64*h, 60.55*m, 299*ns, 3.053*m];
            
            % Build Branching Matrix (Size 12x12)
            n = length(names);
            B = zeros(n-1, n); % Parents x All
            
            % Linear parts
            B(1,2)=1; B(2,3)=1; B(3,4)=1; B(4,5)=1; B(5,6)=1; B(6,7)=1; B(7,8)=1; B(8,9)=1;
            
            % THE BRANCHING: Bi-212 (Idx 9)
            % 64.06% decays to Tl-208 (Idx 11) via Alpha
            % 35.94% decays to Po-212 (Idx 10) via Beta
            B(9, 11) = 0.3594; 
            B(9, 10) = 0.6406;
            
            % Both daughters go to Pb-208 (Idx 12)
            B(10, 12) = 1; % Po-212 -> Pb-208
            B(11, 12) = 1; % Tl-208 -> Pb-208
            
        case 2 % NEPTUNIUM SERIES (4n + 1)
            desc = 'Neptunium Series (Np-237)';
            % Note: Ends in Tl-205 or Bi-209 depending on timeframe. Assuming Bi-209 as quasi-stable.
            names = {'Np-237', 'Pa-233', 'U-233', 'Th-229', 'Ra-225', 'Ac-225', ...
                     'Fr-221', 'At-217', 'Bi-213', 'Ti-209', 'Pb-209', 'Bi-209 (Stable)'};
            
            half_lives = [2.144e6*y, 26.967*d, 1.592e5*y, 7370*y, 14.9*d, 10.0*d, ...
                          4.9*m, 32.3*ms, 45.59*m, 2.161*m, 3.253*h];
            
            % Linear approximation for main plot (Ignoring minor branch of Bi-213 < 2%)
            B = []; % Empty means linear default in DecayMatrix
            
        case 3 % URANIUM SERIES (4n + 2)
            desc = 'Uranium Series (U-238)';
            % The most famous one. Major isotopes only for clarity.
            names = {'U-238', 'Th-234', 'Pa-234m', 'U-234', 'Th-230', 'Ra-226', ...
                     'Rn-222', 'Po-218', 'Pb-214', 'Bi-214', 'Po-214', 'Pb-210', ...
                     'Bi-210', 'Po-210', 'Pb-206 (Stable)'};
                 
            half_lives = [4.468e9*y, 24.1*d, 1.17*m, 2.455e5*y, 7.538e4*y, 1600*y, ...
                          3.8235*d, 3.10*m, 26.8*m, 19.9*m, 164.3*us, 22.3*y, ...
                          5.015*d, 138.376*d];
            
            B = []; % Linear dominant path
            
        case 4 % ACTINIUM SERIES (4n + 3)
            desc = 'Actinium Series (U-235)';
            names = {'U-235', 'Th-231', 'Pa-231', 'Ac-227', 'Th-227', 'Ra-223', ...
                     'Rn-219', 'Po-215', 'Pb-211', 'Bi-211', 'Tl-207', 'Pb-207 (Stable)'};
                 
            half_lives = [7.04e8*y, 25.52*h, 3.276e4*y, 21.772*y, 18.68*d, 11.43*d, ...
                          3.96*s, 1.781*ms, 36.1*m, 2.14*m, 4.77*m];
            
            % Bi-211 has a minor branch (99.7% alpha, 0.3% beta).
            % We simulate the dominant linear path for plotting clarity.
            B = []; 
            
        otherwise
            error('Invalid selection');
    end
end
