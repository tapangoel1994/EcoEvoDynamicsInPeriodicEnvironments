%% Function takes initial conditions, parameter sets, stochastic filtration conditions and returns the number of cycles it takes (upto Maxcycles) for the virus to fall below critical density threshold.

%%Date updated: 1/31/2025
%%Author: Tapan Goel

%% Inputs:
% MaxCycles- maximum number of cycles to run
% params - life history traits and simulation parameters
% criticaldensitythreshold -  minimum density below which a species is
% treated as extinct.
% q_R - fraction of resources to passage between cycles (usually 0)
% q_S - fraction of resources to passage between cycles (usually 0)
% q_E - fraction of resources to passage between cycles (usually 0)
% q_I - fraction of resources to passage between cycles (usually 0)
% q_L - vector of length MaxCycles containing the fraction of lysogens
%       being passed from one cycle to the next.
% q_V - vector of length MaxCycles containing the fraction of virions
%       being passed from one cycle to the next.
% x0 -  initial state of the system before starting the first cycle x0 = [R0,S0,0,0,0,0,0,0,0,0,V_a,0]
% options -  optional variable for the ode solver.

%% Output:
%numcycles - number of cycles (upto MaxCycles) that it takes for everything
%containing viruses (E,I,L and V) to fall below the critical density
%threshold.

function numcycles = StochasticRealization(MaxCycles,params,criticaldensitythreshold,q_R, q_S,q_E,q_I,q_L,q_V,x0,options)
    
    y0 = x0;
    numcycles = MaxCycles;

    if (length(q_L) ~= MaxCycles) | (length(q_V) ~= MaxCycles)
        error('Filtration strings not the same size as # cycles');
    end

    for iter = 1:MaxCycles
        [t_vals,y] = ode113(@ODE_RSEILV_2Species,params.t_vals,y0,options,params);
        %timeseries = [timeseries;y(end,:)];
        TransferMatrix = diag([q_R q_S q_E q_E q_I q_I q_L(iter,1) q_L(iter,1) q_V(iter,1) q_V(iter,1)]);
        y0 = [x0(1) x0(2) zeros(1,8)] + y(end,:)*TransferMatrix; 

        if sum( y0(3:2:10)<criticaldensitythreshold ) == 4 %% Nothing containing the virus is above the threshold

            numcycles = iter;
            break;
        end
        
    end
    
end