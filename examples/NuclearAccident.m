%% Nuclear Accident Simulator: Source term
% Radiological toxicity comparison of chains: I-131 vs Cs-137

clear; clc; close all;
projectRoot = fileparts(pwd);
addpath(fullfile(projectRoot,'src'));

% 1. Typical Reactor Parameters

Power_Thermal = 3000; % MWth (thermal Megawatts, typical commercial reactor)
Fission_Rate_Per_MW = 3.1e16; % Fisions/seg per MW (aprox)
F_total = Power_Thermal * Fission_Rate_Per_MW; % Total fissions per second

SecondsToYear = 60 * 60 * 24 * 365.25; 
Years = input('Enter number of years of operation of the reactor: '); 
T_operation = SecondsToYear * Years; % Operation time in seconds

% Fission Products (U-235) yields
gamma_131 = 0.0289; % ~2.9% for 131
gamma_137 = 0.0619; % ~6.2% for 137

% 2. Simulation time
s = 1; m = 60; h = 3600; d = 86400; y = 3.154e7;

t_start = 1 * h;
t_end   = y*input("Enter the number of simulation's years: ");
num_points = 3000;
T = logspace(log10(t_start), log10(t_end), num_points);

% 3. Chain setting 1: Sb-131 to Xe-131
% Sb-131 -> Te-131m/Te-131 -> I-131 -> Xe-131m -> Xe-131
names_131 = {'Sb131', 'Te131m', 'Te131', 'I131', 'Xe131m', 'Xe131'};
halfs_131 = [23.03*m, 33.25*h, 25.0*m, 8.02*d, 11.84*d];

% Branching matrix (NuDat 3.0)
B131 = zeros(6,6);
B131(1,2) = 0.065; B131(1,3) = 0.935;   % Sb -> Te_m / Te
B131(2,3) = 0.259; B131(2,4) = 0.741;   % Te_m -> Te / I
B131(3,4) = 1.00;                       % Te -> I
B131(4,5) = 0.0039; B131(4,6) = 0.9961; % I -> Xe_m / Xe
B131(5,6) = 1.00;                       % Xe_m -> Xe

% N0 for Iodine computation
% We assumme that the short-lived isotopes are created and destructed in equilibrium (Sb, Te, I)
% Under this condition, the rate of production R=yield*fussions is equal to
% the decay rate, \lambda_I*N_I
% Is assummed I-131 is dominant at t=0
% In reality, there would also be Sb and Te, but I-131 represents the main
% immediate risk
lambda_I = log(2) / halfs_131(4);
N0_I131_total = (gamma_131 * F_total) / lambda_I;
N0_vec_131 = zeros(6,1);
N0_vec_131(4) = N0_I131_total; 

% 4. Chain 2 setting: 131
% Cs-137 -> Ba-137m -> Ba-137
names_137 = {'Cs137', 'Ba137m', 'Ba137'};
halfs_137 = [30.17*y, 2.55*m];

% Branching matrix
B137 = zeros(3,3);
B137(1,2) = 0.944;  % Cs -> Ba_m 
B137(1,3) = 0.056;  % Cs -> Ba   
B137(2,3) = 1.00;   % Ba_m -> Ba

% N0 for Cesio calculation
% Cesium is not in equilibrium: (T_half = 30 years >> T_op = 1 año).
% Formula: N(t) = (Production / lambda) * (1 - exp(-lambda * t_op))
lambda_Cs = log(2) / halfs_137(1);
% Solution for production rate constant:
N0_Cs137_total = (gamma_137 * F_total / lambda_Cs) * (1 - exp(-lambda_Cs * T_operation));

N0_vec_137 = zeros(3,1);
N0_vec_137(1) = N0_Cs137_total;

% 5. Simulation
fprintf('Simulating I-131 chain (N0 = %.2e átomos)...\n', N0_I131_total);
[~, A_131, ~] = DecayMatrix(N0_vec_131, T, halfs_131,'Nuclei', names_131, 'Branching', B131, 'Plotting', false);

fprintf('Simulating Cs-137 chain (N0 = %.2e átomos)...\n', N0_Cs137_total);
[~, A_137, ~] = DecayMatrix(N0_vec_137, T, halfs_137,'Nuclei', names_137, 'Branching', B137, 'Plotting', false);

fprintf('Done.')

% Get activity
Act_I = A_131(4, :); %I-131
% Beta particles of Cs only survive for a few meters, actual danger lies on Ba
Act_Cs = A_137(1, :) + A_137(2, :); % Cs-137 and Ba-137 are in equilibrium

% 6. Visualization (Log-Log)
Min_Visible = 1e-8; 
Act_I(Act_I < Min_Visible) = NaN;
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% Time in days for X axis
T_days = T / d;

loglog(T_days, Act_I, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Iodo-131 (Thyroid)');
hold on;
loglog(T_days, Act_Cs, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Cesio-137 + Ba-137m (Ground)');

% Analysis

% 1. Emergency zone (10 half lives)
xregion_end = 10*halfs_131(4)/d; 
fill([T_days(1) xregion_end xregion_end T_days(1)], ...
     [1e-8 1e-8 1e12 1e12], 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'none', ...
     'DisplayName', 'Acute Phase (KI pills)');

% 2. Crossover Point
diff_A = abs(Act_I - Act_Cs);
[~, idx_cross] = min(diff_A(T_days < 365)); % In the first year
t_cross = T_days(idx_cross);
y_cross = Act_I(idx_cross);

plot(t_cross, y_cross, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'HandleVisibility', 'off');
text(t_cross, y_cross*0.1, sprintf(' Crossover Point\n \\approx %.1f days', t_cross), ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'Color', 'k');

% 3. Historical benchmark
Conv_Ci = 3.7e10;
L_Chernobyl = 1760e15/Conv_Ci; % PBq to Curies, from unscear
% https://www.unscear.org/unscear/uploads/documents/publications/UNSCEAR_2000_Annex-J.pdf
L_safe_Cs  = 1e7 / Conv_Ci; % Limit from the AIEA
L_safe_Ba  = 1e6 / Conv_Ci; % Limit from the AIEA
L_safe_I = L_safe_Ba;
L_safe_CsplusBa = L_safe_Cs + L_safe_Ba;
% https://www-pub.iaea.org/MTCD/Publications/PDF/Pub1578_web-57265295.pdf
yline(4.76e7, '--k', 'Chernobyl', 'LabelHorizontalAlignment', 'right','HandleVisibility', 'off');
yline(1.3e7, '--k', 'Fukushima', 'LabelHorizontalAlignment', 'right','HandleVisibility', 'off');
yline(L_safe_I, '--k', 'Exemption Level of I-131 (IAEA)', 'LabelHorizontalAlignment', 'left','HandleVisibility', 'off');
yline(L_safe_CsplusBa, '--k', 'Excention Level of Cs-137+Ba137m (IAEA)', 'LabelHorizontalAlignment', 'left','HandleVisibility', 'off');

% tags
grid on;
xlabel('Elapsed Time', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Activity (Curies)', 'FontSize', 12, 'FontWeight', 'bold');
title('Source Term Time Evolution: I-131 vs Cs-137', 'FontSize', 14);
subtitle(sprintf('Simulation based on %d years of operations at %d MWth', Years ,Power_Thermal));

legend('Location', 'southwest', 'FontSize', 10);
xlim([T_days(1), T_days(end)]);
ylim([L_safe_I, max([Act_I, Act_Cs])*5]);

% Custom Ticks to understand scale
xticks([1, 7, 30, 90, 365, 365*10, 365*50, 365*350, 356*1050]); % Days
xticklabels({'1 d', '1 w', '1 m', '3 m', ...
    '1 y', '10 y', '50 y','350 y','1050 y'});

y_ticks = 10.^(- 6: 2 : 8); 
yticks(y_ticks);
ytickformat('%.0e');
hold off;