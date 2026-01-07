function [N, A, Lambda] = DecayMonteCarlo(N0, T, half_lives, varargin)
%DECAYMONTECARLO Stochastic simulation of decay chains (Monte Carlo)
%
%   [N, A, Lambda] = DecayMonteCarlo(N0, T, half_lives, Name, Value)
%
%   Simulates the random nature of radioactive decay using Binomial/Multinomial
%   distributions. It adapts to the time steps provided in T (Linear or Log).
%
%   INPUTS:
%       N0          : (Scalar or Vector) Initial number of atoms for the first nucleus.
%       T           : (Vector) Time steps.
%       half_lives  : (Vector) Half lifes.
%
%   OPTIONAL INPUTS:
%       'Branching' : (Matrix) Topology.
%       'Nuclei'    : (Cell) Names for the legend.
%       'TimeUnits' : (Char) Time units (s, m, h, y).
%       'Plotting'  : (Logical) Show plot (default true).
%
%   OUTPUTS:
%       N      : Population [num_nuclei x time] (Integer values)
%       A      : Deterministic Activity (Ci) based on current N.
%       Lambda : Rate matrix (useful for checking topology).

%% 1. Input Parsing
p = inputParser;
addRequired(p, 'N0');
addRequired(p, 'T');
addRequired(p, 'half_lives', @isnumeric);
addParameter(p, 'Branching', [], @isnumeric);
addParameter(p, 'Nuclei', {});
addParameter(p, 'TimeUnits', 's');
addParameter(p, 'Plotting', true, @islogical);

parse(p, N0, T, half_lives, varargin{:});

branch_matrix = p.Results.Branching;
names = p.Results.Nuclei;
units = p.Results.TimeUnits;
plotting = p.Results.Plotting;

%% 2. Setup Physics Constants
num_nuclei = length(half_lives) + 1;
lambdas = log(2) ./ half_lives;

% We build Lambda matrix just for consistency in output/topology check
Lambda = zeros(num_nuclei, num_nuclei);
for i = 1:num_nuclei-1
    Lambda(i, i) = -lambdas(i);
    if isempty(branch_matrix)
        Lambda(i+1, i) = lambdas(i); 
    else
        % Fill off-diagonal for branching 
        daughters = find(branch_matrix(i, :) > 0);
        for d = daughters
            Lambda(d, i) = lambdas(i) * branch_matrix(i, d);
        end
    end
end

%% 3. Stochastic Simulation Loop
steps = length(T);
N = zeros(num_nuclei, steps);

% Initial Condition
current_N = zeros(num_nuclei, 1);
if isscalar(N0)
    current_N(1) = floor(N0);                                               % Ensure integer for Monte Carlo
else
    current_N = floor(N0(:));
end
N(:, 1) = current_N;

% Determine if T is uniform (for plotting logic later)
if steps > 1
    dt_vec = diff(T);
    is_uniform = all(abs(dt_vec - dt_vec(1)) < 1e-9 * dt_vec(1));
else
    is_uniform = true;
end

% --- Main Loop ---
for k = 1:steps-1
    % Calculate time step for this specific interval
    % (Crucial for Logarithmic scales where dt changes)
    dt = T(k+1) - T(k);
    
    if dt <= 0
        N(:, k+1) = current_N;
        continue; 
    end
    
    % Decay probability for this interval: P = 1 - exp(-lambda * dt)
    probs_decay = 1 - exp(-lambdas * dt);
    
    decays = zeros(num_nuclei, 1);
    gains = zeros(num_nuclei, 1);
    
    for i = 1:num_nuclei-1
        if current_N(i) > 0
            % 1. How many decay? (Binomial Distribution)
            % Using MATLAB's Statistics Toolbox function
            num_decayed = binornd(current_N(i), probs_decay(i));
            
            decays(i) = num_decayed;
            
            % 2. Where do they go?
            if num_decayed > 0
                if isempty(branch_matrix)
                    % Linear Chain
                    gains(i+1) = gains(i+1) + num_decayed;
                else
                    % Branching Chain
                    daughters = find(branch_matrix(i, :) > 0);
                    if isempty(daughters)
                        % Default to next if defined in input but not in matrix
                        if i+1 <= num_nuclei
                            gains(i+1) = gains(i+1) + num_decayed;
                        end
                    elseif length(daughters) == 1
                        gains(daughters) = gains(daughters) + num_decayed;
                    else
                        % Multinomial Distribution for branching
                        ratios = branch_matrix(i, daughters);
                        prob_branches = ratios / sum(ratios); % Normalize
                        
                        % mnrnd returns the count for each category
                        branches_counts = mnrnd(num_decayed, prob_branches);
                        
                        for d_idx = 1:length(daughters)
                            target = daughters(d_idx);
                            gains(target) = gains(target) + branches_counts(d_idx);
                        end
                    end
                end
            end
        end
    end
    
    % Update State
    current_N = current_N - decays + gains;
    N(:, k+1) = current_N;
end

%% 4. Compute Activity
% Unit conversion
switch lower(units)
    case {'s', 'sec'}, tf = 1;
    case {'m', 'minute'}, tf = 60;
    case {'h', 'hour'}, tf = 3600;
    case {'y', 'year'}, tf = 3.154e7;
    otherwise, tf = 1; warning('Unknown unit, assuming seconds');
end

lambdas_full = [lambdas(:); 0]; 
% A = lambda * N (in Curies). Note: N is integer here, so A will be "noisy"
A = (lambdas_full ./ tf) .* N ./ 3.7e10;

%% 5. Autoadaptative Visualization (Identical to DecayMatrix)
if plotting && ~isempty(names)
    figure('Color','w', 'Name', 'Monte Carlo Simulation');
    
    % --- Selection of the plotting ---
    if is_uniform 
        plot_cmd = @stairs; % 'stairs' looks better for discrete MC steps
        x_label_txt = ['Time (' units ')'];
        x_limits = [T(1) T(end)];
    else 
        % For log scale, stairs can look weird, standard plot/semilogx is often clearer
        plot_cmd = @semilogx; 
        x_label_txt = ['Time (' units ') [Log Scale]'];
        
        if T(1) == 0
            x_min = T(2); 
        else
            x_min = T(1);
        end
        x_limits = [x_min T(end)];
    end
    
    % Subplot 1: Population
    ax1 = subplot(2,1,1);
    plot_cmd(T, N, 'LineWidth', 1.5); 
    grid on; 
    title('Stochastic Population N(t)');
    ylabel('Nuclei (Count)'); 
    xlim(x_limits);
    hold(ax1,'on')
    
    % VISUAL MARKERS (T_1/2) 
    for i = 1:length(half_lives)
        t_half = half_lives(i);      
        if t_half >= x_limits(1) && t_half <= x_limits(2)
            if ~isempty(names) && i <= length(names)
                lbl = sprintf('T_{1/2}(%s)', names{i});
            else
                lbl = sprintf('T_{1/2}(\\lambda_{%d})', i);
            end   
            xline(ax1, t_half, '--', lbl, ...
                'Color', [0.4 0.4 0.4], ...
                'LabelVerticalAlignment', 'bottom', ...
                'LabelOrientation', 'horizontal', ... 
                'FontSize', 8, ...
                'HandleVisibility', 'off');
        end
    end
    
    if ~isempty(names), legend(names, 'Location', 'best'); end
    hold(ax1, 'off')
    
    % Subplot 2: Activity
    subplot(2,1,2);
    plot_cmd(T, A, 'LineWidth', 1.5);
    grid on; 
    title('Stochastic Activity A(t)');
    ylabel('Activity (Ci)'); 
    xlabel(x_label_txt);
    xlim(x_limits);
end
end