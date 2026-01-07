function [sol] = DecayChain(N_0,T,meanlifes, varargin)
%The function **DecayChain** computes the number of nuclei as a function of time for every nucleus in a
%decay chain, plotting the graphs in a single figure and giving a cell
%'sol' containing these solutions to evaluate at any required time. 
%The chain could has branching or not. 
%To perform the solutions you may need a different number of inputs. 
%
%No branching case
%
%N_0: Number of atoms present at t=0 of the first nucleus of the chain.
%T: Vector containing the initial time, increment and final time for the
%graphs. 
%meanlifes: Vector containing the meanlifes for every nucleus in the right
%order. 
%
%You could provide 2 adittional inputs: nucleos and timeunits (in that orden).
%nucleos: Cell array with elements class char tagging the nuclei. 
%timeunits: Char input indicating in parenthesis the time units in which the meanlifes are given.
%
%Branching case
%
%It works the same way, but the order for the meanlilfes is assignated
%indicating parent, alpha decay product, beta decay product and grandaghter. 
%If it's a gamma decay, you indicate firts the excited state.
%You'll also need to provide another 4 inputs in the next order, but before
%nucleos and timeunits.
%
%Bp: Cell array with vector elements indicating the position of each
%daughter [parent, alpha product, beta product] or similarly for gamma
%decays. 
%Br: Vector ordering the ratio of alpha products by beta products, or similarly for gamma
%decays. 
%Dp: Vector indicating in order the position of the daughter with two
%parents.
%Fp: Cell array with vector elements containing both parents for each
%corresponding daughter in the order [alpha product parent, beta product
%parent] or similarly for gamma decays. 

%% Possibility of linear chains

inputs = length(varargin); % It tells us how many optional inputs we have 

if inputs == 0 || inputs == 1 || inputs == 2 % It's a linear chain without/with nuclei tags input
    linear_chain = true;
    if inputs == 0 % No nuclei tags
        theres_nuclei = false;
        theres_timeunits = false;
    elseif inputs == 1 %There could be the nuclei tags, or the time units of the meanlifes
        if class(varargin{1}) == 'cell' %#ok<BDSCA> % There's only the tags
            theres_nuclei = true;
            theres_timeunits = false;
            nucleos = varargin{1};
        elseif class(varargin{1}) == 'char' %#ok<BDSCA> %There's only the time units
            theres_nuclei = false;
            theres_timeunits = true;
            timeunits = varargin{1};
        end
    else % There are both: nuclei tags and time units
        theres_nuclei = true;
        theres_timeunits = true;
        nucleos = varargin{1};
        timeunits = varargin{2};
    end
elseif 6 - inputs >= 0 % It could be there are not "nucleos"
    if 6 - inputs <= 2 % There has to be at least 4 optionals inputs to compute branching
        linear_chain = false;

        Bp = varargin{1};
        Br = varargin{2};
        Dp = varargin{3};
        Fp = varargin{4};
   
        % Assignation or creation of "nucleos"
        if inputs == 5 %It could be the nuclei tags or time units
            if class(varargin{5}) == 'cell' %#ok<BDSCA> % There's only the tags
                theres_nuclei = true;
                theres_timeunits = false;
                nucleos = varargin{5};
            elseif class(varargin{5}) == 'char' %#ok<BDSCA> %There's only the time units
                theres_nuclei = false;
                theres_timeunits = true;
                timeunits = varargin{5};
            end
        elseif inputs == 6 %There are the both
            theres_nuclei = true;
            theres_timeunits = true;
            nucleos = varargin{5};
            timeunits = varargin{6};
        else % There aren't both: nor nuclei tags nor time units
            theres_nuclei = false;
            theres_timeunits = false;
        end
    else % We have less inputs than needed
        error('To compute a chain with branching there should be given at least Br, Bp, Dp and Fp.')
    end
else % For some reason we have more inputs 
    error('Additional non-required inputs.')
end

%% Preallocation
nuclei = length(meanlifes) + 1 ; % We're addingg + 1 for the missing meanlife of stable element.

y = zeros(1,nuclei); % Decay constants for the nuclei
if ~linear_chain
    no_branchings = length(Bp); % Number of branchings
    y_B = zeros(no_branchings,2); %Dual decay constants y = y_a + y_b for every branching
end

C = zeros(nuclei); % Matrix of recursive constants
N = zeros(nuclei,length(T)); % Matrix of solutions
A = zeros(nuclei,length(T)); % Matrix of acitivities

%% Compute of decay constants
counter = 1;

for ii = 1:nuclei-1
    if linear_chain
        y(ii) = log(2)/meanlifes(ii);
    else
        if counter <= no_branchings
            if ii == Bp{counter}(1)
                y(ii) = log(2)/meanlifes(ii);
                y_B(counter,1) = y(ii) * Br(counter)/(Br(counter)+1); % Alpha decay constant
                y_B(counter,2) = y(ii) - y_B(counter,1); % Beta decay constant
                counter = counter + 1;
            else
                y(ii) = log(2)/meanlifes(ii);
            end
        else
            y(ii) = log(2)/meanlifes(ii);
        end
    end
end

%% Solutions

% Branching parameters
Nbranching = 1;
daughter_counter = 1;

% Creating figure for the graphs
figure(1)
hold on
grid on
if theres_timeunits
    space = ' ';
    time = 'Tiempo';
    xtag = append(time,space,timeunits);
    
    xlabel(xtag,FontName='Times', FontSize=14);
else
    xlabel('Tiempo',FontName='Times', FontSize=14)
end
xlim([0 T(end)])
ylabel('Núcleos',FontName='Times', FontSize=14)
ylim([0 N_0])
%ylim([0 03e15])
legends = cell(nuclei,1);

% Compute of solutions
sol = cell(nuclei,2); %Prealloc to evaluate the solution at any t
for ii = 1:nuclei % Running every nuclei
    solution = zeros(ii,length(T)); %Every row corresponds to the prealloc for partial solutions of N(ii) as a function of time
    for jj = 1:ii % Running only needed constants
        if linear_chain % It's the simple chain without branching
            if ii == 1
                C(ii,jj) = N_0; %Initial condition for the parent
            else
                if ii - jj > 0
                    C(ii,jj) = C(ii-1,jj) * y(ii-1)/(y(ii) - y(jj)); %Recursive relation
                else
                    C(ii,jj) = -sum(C(ii,:)); %Initial condition
                end
            end
        else %The chain has branching
            if ii == 1 %
                C(ii,jj) = N_0; % 
            elseif Nbranching <= length(Bp) % There are branchings left in the chain
                if ii == Dp(daughter_counter) % It's a daughter with two parents
                    ii_f1 = Fp{daughter_counter}(1); % Beta parent with alpha decay
                    ii_f2 = Fp{daughter_counter}(2); % Alpha parent with beta decay
                    if ii_f1 == Bp{Nbranching}(1) % Parent from alpha decay presents branching
                        if ii - jj > 0
                            C(ii,jj) = C(ii_f1,jj) * y_B(Nbranching,2)/(y(ii) - y(jj)) + C(ii_f2,jj) * y(ii_f2)/(y(ii) - y(jj));
                        else
                            C(ii,jj) = -sum(C(ii,:));
                            daughter_counter = daughter_counter + 1;
                        end
                    elseif ii_f2 == Bp{Nbranching}(1) % Parent form beta decay presents branching
                        if ii - jj > 0
                            C(ii,jj) = C(ii_f1,jj) * y(ii_f1)/(y(ii) - y(jj)) + C(ii_f2,jj) * y_B(Nbranching,1)/(y(ii) - y(jj));
                        else
                            C(ii,jj) = -sum(C(ii,:));
                            daughter_counter = daughter_counter + 1;
                        end
                    else % Neither of both parents has branching
                        if ii - jj > 0
                            C(ii,jj) = C(ii_f1,jj) * y(ii_f1)/(y(ii) - y(jj)) + C(ii_f2,jj) * y(ii_f2)/(y(ii) - y(jj));
                        else
                            C(ii,jj) = -sum(C(ii,:));
                            daughter_counter = daughter_counter + 1;
                        end
                    end
                elseif ii == Bp{Nbranching}(2) % Daughter by alpha emission
                    ii_f = Bp{Nbranching}(1); % Alpha emitter parent
                    if ii - jj > 0
                        C(ii,jj) = C(ii_f,jj) * y_B(Nbranching,1)/(y(ii) - y(jj));
                    else
                        C(ii,jj) = -sum(C(ii,:));
                        if ii+1 == Dp(daughter_counter)
                            Nbranching = Nbranching + 1; % We must increase Nbranching if the next nuclei has two parents
                        end
                    end
                elseif ii == Bp{Nbranching}(3) % Daughter by beta emission
                    ii_f = Bp{Nbranching}(1); % Beta emitter parent
                    if ii - jj > 0 
                        C(ii,jj) = C(ii_f,jj) * y_B(Nbranching,2)/(y(ii) - y(jj));
                    else
                        C(ii,jj) = -sum(C(ii,:));
                        Nbranching = Nbranching + 1; % We must increase the Nbranching
                    end
                else % This nuclei has only 1 parent, but there are others ahead with branching
                    if ii - jj > 0
                        C(ii,jj) = C(ii-1,jj) * y(ii-1)/(y(ii) - y(jj));
                    else
                        C(ii,jj) = -sum(C(ii,:));
                    end
                end                
            else % Linear chain again or gamma decay in the first case.
                if ii == Dp(end) % It's the last daughter with two parents
                    ii_f1 = Fp{length(Dp)}(1);
                    ii_f2 = Fp{length(Dp)}(2);
                    if ii_f1 == Bp{Nbranching-1}(1) % Gamma decay in the first case
                        if ii - jj > 0
                            C(ii,jj) = C(ii_f1) * y_B(Nbranching-1,2)/(y(ii) - y(jj)) + C(ii_f2,jj) * y(ii_f2)/(y(ii) - y(jj));
                        else
                            C(ii,jj) = -sum(C(ii,:));
                        end
                    else % The nuclei has 2 parents decaying linearly
                        if ii - jj > 0
                            C(ii,jj) = C(ii_f1,jj) * y(ii_f1)/(y(ii) - y(jj)) + C(ii_f2,jj) * y(ii_f2)/(y(ii) - y(jj));
                        else
                            C(ii,jj) = -sum(C(ii,:));
                        end
                    end
                else % The nuclei has only 1 parent
                    if ii - jj > 0
                        C(ii,jj) = C(ii-1,jj) * y(ii-1)/(y(ii) - y(jj));
                    else
                        C(ii,jj) = -sum(C(ii,:));
                    end
                end
            end
        end
        eq = @(t) C(ii,jj) * exp(-y(jj)*t); % Equation for the constant C(ii,jj)
        solution(jj,:) = eq(T); % Partial solution as function of time (Only the apport of this constant to the whole solution)
    end

    %% Graphs

    %Molibdeno
    %mm = 99;
    if ii == 1
        N(ii,:) = solution;
        A(ii,:) = (1/((3600)*(3.7e10))) * y(ii) * N(ii,:); %Actividad en curies
        plot(T,N(ii,:),LineWidth=2)
        if theres_nuclei
            legends{ii} = sprintf(nucleos{ii});
        else
            legends{ii} = sprintf('\\lambda_{%d}',ii);
        end
        sol{ii,1} = @(t) sum(C(ii,:) .* exp(-y*t)); %N(ii) as a function of time.
        sol{ii,2} = @(t) (1/((3600)*(3.7e10))) * y(ii) * sum(C(ii,:) .* exp(-y*t)); % A(ii) as a function of time in seconds
        %We have the activity in Curies (Ci), i.e. desintegrations per s when we multiply the activity
        %in desintegrations per hour by the factor (1/(3600 * 3.7x10^10))
    else
        N(ii,:) = sum(solution);
        A(ii,:) = (1/((3600)*(3.7e10))) * y(ii) * N(ii,:);
        plot(T,N(ii,:),LineWidth=2)
        if theres_nuclei
            legends{ii} = sprintf(nucleos{ii});
        else
            legends{ii} = sprintf('\\lambda_{%d}',ii);
        end
        sol{ii,1} =@(t) sum(C(ii,:) .* exp(-y*t));
        sol{ii,2} = @(t) (1/((3600)*(3.7e10))) * y(ii) * sum(C(ii,:) .* exp(-y*t));
    end
end

legend(legends);
hold off

%% Activity Graphs

figure(2)
hold on 
grid on
if theres_timeunits
    space = ' ';
    time = 'Tiempo';
    xtag = append(time,space,timeunits);
    
    xlabel(xtag,FontName='Times', FontSize=14);
else
    xlabel('Tiempo',FontName='Times', FontSize=14)
end
xlim([0 T(end)])
ylabel('Actividad (Ci)',FontName='Times', FontSize=14)
ylim([0 max(max(A))])

for ii = 1:nuclei
    plot(T,A(ii,:),LineWidth=2)
end

legend(legends)

end