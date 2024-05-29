%% Code implements single time step for a RSEILV model with no MOI dependence.  
% Modified from ODE_RSEILV_2species
%% 
%%Author: Tapan Goel
%%Date of creation: 05/29/2024

%%Input:
% t = timepoint
% y = (R(t),S(t),E1(t),E2(t),I1(t),I2(t),L1(t),L2(t),V1(t),V2(t)) : system CONCENTRATION state at timepoint t.
% params = structure of parameter values with.

%%Output:
% dydt = d(R(t),S(t),E1(t),E2(t),I1(t),I2(t),L1(t),L2(t),V1(t),V2(t))/dt : the derivative of the system state at timepoint t

function dydt = ODE_RSEILV_2species_ResourceDependentLysis(t,y,params)

y = reshape(y,[1,10]);
R = y(1);
S = y(2);%*heaviside(y(2)-1/flask_volume);
E = y(3:4);%*heaviside(y(3)-1/flask_volume);
I = y(5:6);%*heaviside(y(4)-1/flask_volume);
L = y(7:8);%*heaviside(y(5)-1/flask_volume);
V = y(9:10);%*heaviside(y(6)-1/flask_volume);

q = reshape(params.q,[1,2]);

psi = params.mu_max*R/(params.R_in+R);
eta = params.eta*R/(params.R_eta + R);

dydt = zeros(1,10);

dydt(1) = - params.conversion_efficiency*psi*(S+sum(E)+sum(I)+sum(L)) - params.d_R*R; %% dR/dt

dydt(2) = psi*S - params.phi*S*sum(V) - params.d_S*S; %% dS/dt

dydt(3:4) = params.phi*S.*V - params.lambda*E - params.d_E*E; %% dE/dt

dydt(5:6) = (1-reshape(params.q,[1,2])).*params.lambda.*E - eta*I + reshape(params.gamma,[1,2]).*L - params.d_I*I; %% dI/dt

dydt(7:8) = psi*L + reshape(params.q,[1,2]).*params.lambda.*E - reshape(params.gamma,[1,2]).*L - params.d_L*L; %% dL/dt

dydt(9:10) = params.bet*eta*I - params.phi*V*( S+sum(E+I+L) ) - params.m*V; %% dV/dt 


dydt = reshape(dydt,[10,1]);
        
end

