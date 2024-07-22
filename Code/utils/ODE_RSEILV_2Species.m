%% Code implements single time step for a RSEILV model with no MOI dependence corresponding to eq 7 in the main manuscript

%%Author: Tapan Goel
%%Date of creation: 9/14/2023

%%Input:
% t = timepoint
% y = (R(t),S(t),Ea(t),Eb(t),Ia(t),Ib(t),La(t),Lb(t),Va(t),Vb(t)) : system CONCENTRATION state at timepoint t.
% params = structure of parameter values

%%Output:
% dydt = d(R(t),S(t),Ea(t),Eb(t),Ia(t),Ib(t),La(t),Lb(t),Va(t),Vb(t))/dt : the derivative of the system state at timepoint t

function dydt = ODE_RSEILV_2species(t,y,params)

y = reshape(y,[1,10]);
R = y(1);
S = y(2);%*heaviside(y(2)-1/flask_volume);
E = y(3:4);%*heaviside(y(3)-1/flask_volume);
I = y(5:6);%*heaviside(y(4)-1/flask_volume);
L = y(7:8);%*heaviside(y(5)-1/flask_volume);
V = y(9:10);%*heaviside(y(6)-1/flask_volume);

p = reshape(params.p,[1,2]); %integration probabilities

psi = params.mu_max*R/(params.R_in+R); %monod resource consumption

dydt = zeros(1,10);

dydt(1) = - params.conversion_efficiency*psi*(S+sum(E)+sum(I)+sum(L)) - params.d_R*R; %% dR/dt

dydt(2) = psi*S - params.phi*S*sum(V) - params.d_S*S; %% dS/dt

dydt(3:4) = params.phi*S.*V - params.lambda*E - params.d_E*E; %% dE/dt

dydt(5:6) = (1-reshape(params.p,[1,2])).*params.lambda.*E - params.eta*I + reshape(params.gamma,[1,2]).*L - params.d_I*I; %% dI/dt

dydt(7:8) = psi*L + reshape(params.p,[1,2]).*params.lambda.*E - reshape(params.gamma,[1,2]).*L - params.d_L*L; %% dL/dt

dydt(9:10) = params.bet*params.eta*I - params.phi*V*( S+sum(E+I+L) ) - params.m*V; %% dV/dt 


dydt = reshape(dydt,[10,1]);
        
end

