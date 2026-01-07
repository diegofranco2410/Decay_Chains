function [N, A, Lambda] = DecayMatrix(N0, T, half_lives, varargin)
%DECAYMATRIX solves decay chains using exponential matrix
%
%   [N, A, Lambda] = DecayMatrix(N0, T, half_lives, Name, Value)
%
%   The functions buils the transition rates matrix (lambda) and solves the
%   system of differential equations dN/dt=lambda * N using expm()
%
%   INPUTS:
%       N0          : (Scalar or Vector) Initial quantity. If scalar, it's assumed
%                     it represents the first nucleus.
%       T           : (Vector) Time.
%       half_lives  : (Vector) Half lifes.
%
%   OPTIONAL INPUTS:
%       'Branching' : (Matrix) Conexion matrix (Topology).
%                     Branching(i, j) = Probability that i decays into j. 
%                     If it's not given, a linear chain is assumed: 1->2->3...
%       'Nuclei'   : (Cell) Names for the legend.
%       'TimeUnits' : (Char) Time units for half-lifes: 
%                     seconds (s), minutes (m), hours (h) or years(y)
%                     If it's not given, seconds are assumed.
%       'Plotting' : Logical value (true or false) specifying if the user
%                    wants to display the plots. If it's not given, true is
%                    assumed.
%
%   OUTPUTS:
%       N      : Populatio [num_nuclei x time]
%       A      : Activity (Ci) [num_nuclei x time]
%       Lambda : Rate matrix.

%% 1. Input Parsing
p = inputParser;
addRequired(p, 'N0');
addRequired(p, 'T');
addRequired(p, 'half_lives', @isnumeric);                                   % Check if the input is numeric
addParameter(p, 'Branching', [], @isnumeric);                               % If not given, assign the empty matrix, and check the input is numeric
addParameter(p, 'Nuclei', {});                                              % If not give, asign the empty cell for legends
addParameter(p, 'TimeUnits', 's');                                          % If not given, asign 's' (secods) as the unit
addParameter(p, 'Plotting', true, @islogical);                              % By default, it shows the plot

parse(p, N0, T, half_lives, varargin{:});                                   % Check inputs 

branch_matrix = p.Results.Branching;
names = p.Results.Nuclei;
units = p.Results.TimeUnits;
plotting = p.Results.Plotting;

%% 2. Construction of Lambda Matrix 

n = length(half_lives) + 1;                                             % Need to add the stable one

lambdas = log(2) ./ half_lives;                                             % Decay constants vector
Lambda = zeros(n, n);                                                       % Prealloc 

% Rule: Lambda(row, col) -> rate row "col" turning into "row".
%        Lambda(col, col)  -> total vanishing rate of col.

% a) Fill the main diagonal (vanish)
for i = 1:n-1
    Lambda(i, i) = -lambdas(i);
end
% Last element is stable, Lambda(n,n) = 0.

% b) Fill outside of the main diagonal (Production)
if isempty(branch_matrix)
    % Default case: linear chain (1->2, 2->3, ...)
    for i = 1:n-1
        Lambda(i+1, i) = lambdas(i);                                        % Everything decaying from 'i' goes to 'i+1'
    end
else
    % Branching case: We use the branching matrix from the user 
    % branch_matrix(i, j) is the fraction of i that goes to j.
    % Note: branch_matrix is squared of the size of parents
    [rows, cols] = size(branch_matrix);
    if rows > n || cols > n
        error('The branching matrix if larger than the number of parents.');
    end
    
    for i = 1:rows       % Parent (Source)
        for j = 1:cols   % Daughter (Target)
            fraction = branch_matrix(i, j);
            if fraction > 0
                % Production rate of the daughter 'j' coming from father 'i'
                % is: lambda_parent * fraction.
                % In the lambda matrix, target is row, source is col.
                Lambda(j, i) = lambdas(i) * fraction;
            end
        end
    end
end

%% 3. Solution with Exponential Matrix
% Initial Conditions Vector
N0_vec = zeros(n, 1); 
if isscalar(N0)                                                             % Only the initial parent has initial nuclei
    N0_vec(1) = N0;                                                         % First entry N0, all others 0
else
    N0_vec = N0(:);                                                         % Copy of the N0 vector of the user
end

N = zeros(n, length(T));                                                    % Prealloc for N(t)

% Determine if the step in time is constant (for future optimization)
if length(T) > 1
    dt_vec = diff(T);
    dt = dt_vec(1);
    is_uniform = all(dt_vec - dt == 0);
else % Trivial case, only one point
    is_uniform = true;
    dt = 0;
end

% Temporal Loop
if is_uniform && length(T)>1 % Optimizd method
    % N(t+dt) = EN(t)
    E = expm(Lambda * dt);
    if T(1) == 0 
        N(:,1) = N0_vec;                                                    % Set initial population at time t=0
    else
        N(:,1) = expm(Lambda * T(1)) * N0_vec;
    end

    for k = 2:length(T)
        N(:, k) = E * N(:, k-1);                                            % Update population for the next time step
    end
else % Robust method (for logarithmic scales or arbitrary time)
    for t_idx = 1:length(T)
        t = T(t_idx);
        N(:,t_idx) = expm(Lambda * t) * N0_vec;
    end
end

%% 4. Compute Activity
% Unit conversion to Curies 
switch lower(units)
    case {'s', 'sec'}, tf = 1;
    case {'m', 'minute'}, tf = 60;
    case {'h', 'hour'}, tf = 3600;
    case {'y', 'year'}, tf = 3.154e7;
    otherwise, tf = 1; warning('Unidad desconocida, asumiendo segundos');
end

% A = lambda * N (Converted to Ci)
% Note: vector 'lambdas'.
lambdas_full = [lambdas(:); 0];                                             % Adding the decay constant 0 for last stable nucleus
A = (lambdas_full ./ tf) .* N ./ 3.7e10;

%% 5. Autoadaptative Visualization
if plotting && ~isempty(names)
    figure('Color','w');
    
    % --- Selection of the plotting ---
    if is_uniform % Equally spaced time
        plot_cmd = @plot;                                                   % Pointer to plot
        x_label_txt = ['Time (' units ')'];
        x_limits = [T(1) T(end)];
    else 
        plot_cmd = @semilogx;                                               % Pointer to semilogx
        x_label_txt = ['Time (' units ') [Log Scale]'];
        
        % Adjuts limits to avoid log(0)
        if T(1) == 0
            x_min = T(2); 
        else
            x_min = T(1);
        end
        x_limits = [x_min T(end)];
    end
    
    % Subplot 1: Population
    ax1 = subplot(2,1,1);
    plot_cmd(T, N, 'LineWidth', 2); 
    grid on; 
    title('Populatio N(t)');
    ylabel('Nuclei'); 
    xlim(x_limits);
    hold(ax1,'on')
    % VISUAL MARKERS (T_1/2) 
    for i = 1:length(half_lives)
        t_half = half_lives(i);      
        % Make sure it's within the x axis
        if t_half >= x_limits(1) && t_half <= x_limits(2)
            % Generate tag (Name or index)
            if ~isempty(names) && i <= length(names)
                lbl = sprintf('T_{1/2}(%s)', names{i});
            else
                lbl = sprintf('T_{1/2}(\\lambda_{%d})', i);
            end   
            % Vertical line using xline 
            % HandleVisibility='off' so it doesn't show in legend
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
    
    % Subplot 2: Actividad
    subplot(2,1,2);
    plot_cmd(T, A, 'LineWidth', 2);
    grid on; 
    title('Activity A(t)');
    ylabel('Activity (Ci)'); 
    xlabel(x_label_txt);
    xlim(x_limits);
end

end