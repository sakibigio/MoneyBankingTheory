function [R_f]=averagerate_leontief(theta_0,delta_r,matchtech);
lambda= matchtech.lambda; % Matching Efficienty
eta   = matchtech.eta;

% Inputs: theta_0, lambda, t in [0,1]
    
% Defines the case in which theta_0=1
if theta_0 == 1
    R_f=(1-eta)*delta_r;
end

% Defines the case in which theta_0>1
if theta_0 > 1
    bartheta=1+(theta_0-1)*exp(lambda);
end

% Defines the case in which theta_0<1
if theta_0 < 1
    bartheta=1/(1+(theta_0^(-1)-1)*exp(lambda));
end

R_f=delta_r*(bartheta-bartheta^(eta)*theta_0^(1-eta))/(bartheta-theta_0);


% % Test 
% check =theta_0>((exp(lambda)-1)/(exp(lambda)+1))^2;
% check2=theta_0<((exp(lambda)+1)/(exp(lambda)-1))^2;
% T_ok=check*check2;