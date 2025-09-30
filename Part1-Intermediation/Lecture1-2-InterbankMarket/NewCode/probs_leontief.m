function [Gammad,Gammas]=probs_leontief(theta_0,matchtech);
lambda=matchtech.lambda; % Matching Efficienty

% Inputs: theta_0, lambda, t in [0,1]
    
% Defines the case in which theta_0=1
if theta_0 == 1
    Gammas=1-exp(-lambda);
    Gammad=1-exp(-lambda);
end

% Defines the case in which theta_0>1
if theta_0 > 1
   % theta1  = 1 + (theta_0-1)*exp(lambda);
    Gammas=1-exp(-lambda);
    Gammad=Gammas/theta_0;
end

% Defines the case in which theta_0<1
if theta_0 < 1
    Gammad=1-exp(-lambda);
    Gammas=Gammad*theta_0;
end


% % Test 
% check =theta_0>((exp(lambda)-1)/(exp(lambda)+1))^2;
% check2=theta_0<((exp(lambda)+1)/(exp(lambda)-1))^2;
% T_ok=check*check2;