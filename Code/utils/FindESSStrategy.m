% Function takes PIP as input and outputs the ESS strategy if one exists

%%Date Created: 1/09/2025
%%Date updated: 1/31/2025
%%Author: Tapan Goel

%% Inputs:
% PIP - a square matrix whose (i,j)'th element contains the result of an
%       invasion by the phenotype InvasionVariable(j) in a resident population of
%       phenotype InvasionVariable(i);
% InvasionVariable - vector of values of the variable across which invasion
% dynamics have been ascertained.
%% Output:
% ESS - Evolutionary stable strategy value of the invasion variable. This
%       is a value of the invasion variable such that it can invade any other strategy but is resistant to invasion by other strategies.

function ESS = FindESSStrategy(PIP,InvasionVariable)

    n = length(InvasionVariable);

    ESS = [];

    for strategynum = 1:n
        
       
        if sum(PIP(:,strategynum) == -2) == n  %% check if strategy can be resident in the long term
            %sprintf("Strategy P=%.2f is not long term viable",InvasionVariable(strategynum))

        else %% if strategy is viable in the long term

            not_invadable = 0;
            always_invades = 0;
            %% Check that if strategy num is resident, no one can displace it
        
            if sum(PIP(:,strategynum) == 1) == 0
                not_invadable = 1;
            end

            %% Check that if strategy is mutant, it can displace everyone
            if sum(PIP(strategynum,:) == -1) == 0
                always_invades = 1;
            end

            if not_invadable == 1 & always_invades == 1
                ESS = [ESS; InvasionVariable(strategynum)];
            end
        end
    end
end


               

        

