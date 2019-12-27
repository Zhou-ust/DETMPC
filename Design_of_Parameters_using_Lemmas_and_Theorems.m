% % Application to power system
% This demo shows how to design the appropriate parameters for the power
% systems in the submitted manuscript named: 
% 'Distributed event-triggered predictive control for a class of nonlinear
% interconnected systems'
%
% References
%
% [1] Y. Zhou, D. Li, J. Lu, Y. Xi, and L. Cen. Networked and distributed 
%     predictive control of nonlinear systems subject to asynchronous communication. 
%     IET Control Theory & Applications, vol. 12, no. 4, pp. 504-514, 2017
% 
% [2] Y. Zhou, D. Li, Y. Xi, and Z. Gan. Periodic event-triggered control 
%     for distributed networked multiagents with asynchronous communication: 
%     A predictive control approach. International Journal of Robust and 
%     Nonlinear Control, vol. 29, no. 1, pp. 43-66, 2019.
%
% Copyright 2019 Y. Zhou 
% Contact: zhouyuanq@126.com

clc
clear

% Parameters for the 2-machine power system
omg0 = 50;      %steady state frequency
Nm   = 2;           %# of machines
tour = 0.1;

H = [5   4];
D = [1      1.5];
T = [4      5];
E = [0.5    0.8];
delta0 = [51.0 51.5];
X  = 10;

Qi = [200 0 0;
    0 10 0;
    0 0 50];
Ri = 1;

% Given the dynamics of the 2-machine power system
% The first maching
A1 = [0 -1 0; 0 -D(1)/2/H(1) omg0/2/H(1);  0 0 -1/T(1)];
B1 = [0;0;1/T(1)];
% dlt12 = x1(1) + 51 - x2(1) - 51.5;
% Pe1 = 0.5*0.8/10*(sin(x1(2) - x2(2)) - sin(dlt12));
% d1 = 0; % disturbence

% The second maching
A2 = [0 -1 0; 0 -D(2)/2/H(2) omg0/2/H(2);  0 0 -1/T(2)];
B2 = [0;0;1/T(2)];
% dlt12 = x1(1) + 51 - x2(1) - 51.5;
% Pe1 = 0.5*0.8/10*(sin(x1(2) - x2(2)) - sin(dlt12));
% d2 = 0; % disturbence

% Obtain the kappa1 kappa2
kappa1 = omg0*E(1)*E(2)/X/H(1)
kappa2 = omg0*E(2)*E(1)/X/H(2)

% Design the feedback gains K1 K2
% eigA1 = eig(A1);
% [K1,S1,P1] = lqr(A1,B1,Qi,Ri);
K1 = -lqr(A1,B1,Qi,Ri)
% K1 = [-259.9324 -44.1761 -153.1983]
% K1 = [-31.6228   -34.4080  -112.1878]

Ad1 = A1 + B1*K1;
% eigAd1 = eig(Ad1)

% eigA2 = eig(A2)
K2 = -lqr(A2,B2,Qi,Ri)
% K2 = [-259.9324 -44.1761 -153.1983]
% K2 = [-31.6228   -34.4080  -112.1878]

Ad2 = A2 + B2*K2;
% eigAd2 = eig(Ad2)


% Next, compute the matrix P for Assumption 2.
%  Agent 1;
bar_Q1 = Qi + K1'*Ri*K1;
Lyap1 = Ad1 + kappa1*eye(3);
P1 = lyap(Lyap1', bar_Q1)

%  Agent 2;
bar_Q2 = Qi + K2'*Ri*K2;
Lyap2 = Ad2 + kappa2*eye(3);
P2 = lyap(Lyap2', bar_Q2)

P_test = [P1 zeros(3);
    zeros(3) P2];

% Design the parameters
Alpha1 = 0.9;
Beta1 = 0.7;

Alpha2 = 0.8;
Beta2 = 0.8;

Epsilon1 = 0.03;
Epsilon2 = 0.03;

% Obtains the bounds of the prediction horizon T
% First, using Eq.(16)
T1min  = -2*min(eig(P1))/max(eig(bar_Q1))/Beta1*log(Alpha1);
T2min  = -2*min(eig(P2))/max(eig(bar_Q2))/Beta2*log(Alpha2);
Tmin = max(T1min, T2min)

% Second, test
T= 1

left191 = exp(abs(max(real(eig(Ad1))))*Beta1*T)*Beta1^2*T;
right191 = Alpha1/kappa1;
left192 = exp(max(real(eig(Ad2)))*Beta2*T)*Beta2^2*T;
right192 = Alpha2/kappa2;
if ~isempty(left191) && ~isempty(right191) && ~isempty(left192) && ~isempty(right192)
    if left191 <= right191 && left192 <= right192
    disp(['The conidtion Eq.(16) and Eq. (19) is satisfied.' ]);
    else
        disp(['Choose the prediction horizon T with T >=' num2str(Tmin)]);
    end
end

% End of determining the propriate T 




% Design q using Eq.(20)
q_right201 = Alpha1/(1-Beta1)-1;
q_right202 = Alpha2/(1-Beta2)-1;
q_min = max(q_right201, q_right202);

q = 3
if ~isempty(q_min) && ~isempty(q)
    if q <= q_min
    disp(['The conidtion Eq.(20) is satisfied.' ]);
    else
        disp(['Choose the q with q >=' num2str(q_min)]);
    end
end


% Design rho1 and rho2 using Eq.(8) and Eq.(17)
% For rho1, 
Eta1 = max(eig(A1)) + kappa1*Nm;
rho1max = (q+1 - Alpha1)/max(sqrt(eig(P1)))/(q+1)/Beta1/T *  Epsilon1*...
    exp(- (Eta1*Beta1*T+ max(eig(A1))*(1-Beta1)*T^2)) - ...
    (1 + Alpha1)* Epsilon1 /max(sqrt(eig(P1)))/(q+1) * kappa1*Nm;

rho1 = 0.0004
if ~isempty(rho1max) && ~isempty(rho1)
    if rho1 <= rho1max
    disp(['The conidtion Eq.(17) is satisfied.' ]);
    else
        disp(['Choose the rho1 with rho1 <=' num2str(rho1max)]);
    end
end

% For rho2,
Eta2 = max(eig(A2)) + kappa1*Nm;
rho2max = (q+1 - Alpha2)/max(sqrt(eig(P2)))/(q+1)/Beta2/T *  Epsilon2*...
    exp(- (Eta2*Beta2*T+ max(eig(A2))*(1-Beta2)*T^2)) - ...
    (1 + Alpha2)* Epsilon2 /max(sqrt(eig(P2)))/(q+1) * kappa2*Nm;

rho2 = 0.0001
if ~isempty(rho2max) && ~isempty(rho2)
    if rho2 <= rho2max
    disp(['The conidtion Eq.(17) is satisfied.' ]);
    else
        disp(['Choose the rho2 with rho2 <=' num2str(rho2max)]);
    end
end



% Check the stability condition using Eq.(22)
Sigma1 = (1+Alpha1)* Epsilon1/(q+1)*kappa1*Nm + rho1*max(eig(P1));
Sigma2 = (1+Alpha2)* Epsilon2/(q+1)*kappa2*Nm + rho2*max(eig(P2));
bar_sigma1 = Beta1*Sigma1*T*exp(Eta1*Beta1*T);
bar_sigma2 = Beta2*Sigma2*T*exp(Eta2*Beta2*T);
left221 = min(eig(Qi))/max(eig(P1))*Beta1*T*Beta1*(Epsilon1 - bar_sigma1)^2
w231 = max(sqrt(eig(Qi)))^2/min(sqrt(eig(P1)))^2;
TEMP1 = ((2*bar_sigma1*T*Alpha1*Epsilon1)+(q+1)*bar_sigma1^2)/2/(q+1)/max(abs(eig(A1)));
TEMP2 = exp(2*max(eig(A1))*(1-Beta1)*T)-1;
right221 = w231*(TEMP1* TEMP2 + 2*bar_sigma1*Alpha1*Epsilon1*(1-Beta1)/(q+1)/Beta1)

if ~isempty(left221) && ~isempty(right221)
    if left221 <= right221
    disp(['The stability conidtion Eq.(22) is satisfied.' ]);
    else
        disp(['Choose the parameters again, like T, Q, rho1, and rho2']);
    end
end


% calculate the bounded compact set \Omega(bar_epsilon1) for generator 1; 
bar_epsilon1 = sqrt(2*max(eig(P1))*norm(sqrt(P1),2)/min(eig(Qi))*Epsilon1*rho1)
if ~isempty(bar_epsilon1)
    disp(['The state 1 will converge to the set \Omega(' num2str(bar_epsilon1) ')' ]);
end

% calculate the bounded compact set \Omega(bar_epsilon1) for generator 1; 
bar_epsilon2 = sqrt(2*max(eig(P2))*norm(sqrt(P2),2)/min(eig(Qi))*Epsilon2*rho2)
if ~isempty(bar_epsilon2)
    disp(['The state 2 will converge to the set \Omega(' num2str(bar_epsilon2) ')' ]);
end

























