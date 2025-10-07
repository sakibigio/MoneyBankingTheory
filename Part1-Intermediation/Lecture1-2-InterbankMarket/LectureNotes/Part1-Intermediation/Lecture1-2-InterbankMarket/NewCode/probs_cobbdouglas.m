function [Gammad,Gammas,T_max]=probs_cobbdouglas(theta_0,matchtech);
lambda=matchtech.lambda; % Matching Efficienty

% Inputs: theta_0, lambda, t in [0,1]
sqrt_theta0 = sqrt(theta_0);
alpha = (1 + sqrt_theta0) / (1 - sqrt_theta0);
T_star = (1 / lambda) * log(abs(alpha));
T_max = min(T_star, 1);

Gammas=1-exp(-lambda*T_max)*((alpha+exp(lambda*T_max))/(alpha+1))^2;
Gammad=1-exp(-lambda*T_max)*((alpha-exp(lambda*T_max))/(alpha-1))^2;


% % Test 
% check =theta_0>((exp(lambda)-1)/(exp(lambda)+1))^2;
% check2=theta_0<((exp(lambda)+1)/(exp(lambda)-1))^2;
% T_ok=check*check2;