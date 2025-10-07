function [CHIm_t,CHIp_t,r_f_t,gammad_t,gammas_t,theta_1]=analytic_cobbdouglas(t,theta_0,delta_r,matchtech)
% alpha=matchtech.alpha;   % Technology for matching --- trade breakdown probability
lambda=matchtech.lambda; % Matching Efficienty
eta=matchtech.eta;       % Bargaining Power of Borrower

% Inputs: theta_0, lambda, t in [0,1]
sqrt_theta0 = sqrt(theta_0);
alpha = (1 + sqrt_theta0) / (1 - sqrt_theta0);
T_star = (1 / lambda) * log(abs(alpha));
T_max = min(T_star, 1);

% Test 
check =theta_0>((exp(lambda)-1)/(exp(lambda)+1));
check2=theta_0<((exp(lambda)+1)/(exp(lambda)-1));
T_ok=check*check2;
if theta_0 < 1
    if t <= T_max
        exp_term = alpha * exp(-lambda * t);
        theta_t = ((exp_term - 1) / (exp_term + 1))^2;

    if T_star>1
       exp_T = alpha * exp(-lambda * T_max);
       theta_1 = ((exp_T - 1) / (exp_T + 1))^2;
    else
       theta_1 = 0;
    end
  %     Psi_plus = 1 - ((exp_T + 1) / (alpha + 1))^2;
  %     Psi_minus = 1 - ((exp_T - 1) / (alpha - 1))^2;
    else
        theta_t = 0;
        theta_1 = 0;
        Psi_plus = 0;
        Psi_minus= 1;
    end
elseif theta_0 > 1
    % same formula applies due to symmetry — alpha > 1 in this case too
    if t <= T_max
        exp_term = alpha * exp(-lambda * t);
        theta_t = ((exp_term - 1) / (exp_term + 1))^2;
    if T_star>1
       exp_T = alpha * exp(-lambda * T_max);
       theta_1 = ((exp_T - 1) / (exp_T + 1))^2;
    else
       theta_1 =10e16;
    end

    else
        theta_t = inf;
        theta_1 = inf;
        Psi_plus = 1;
        Psi_minus= 1;
    end
else
    % theta_0 == 1, special case
    theta_t = 1;
    theta_1 = 1;
    Psi_plus  = 1 - exp(-lambda);  % limit as alpha → ∞
    Psi_minus = 1 - exp(-lambda); % same here
end

if theta_0 ~= 1
    CHIp_t = (delta_r)*(theta_1/theta_t)^eta * ((theta_t-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    CHIm_t = (delta_r)*(theta_1/theta_t)^eta * ((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    Sigma_t= CHIm_t-CHIp_t;
    r_f_t  = (theta_1/theta_t)^eta * (((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1))-eta*((1-theta_t)/(1-theta_1)))*(delta_r) ;
else
   Sigma_t= (delta_r)*exp(-lambda * (1-t));
   CHIp_t = (delta_r)*((1-eta)*(1-exp(-lambda * (1-t))));
   % CHIp_t = (delta_r)*((1-eta)*(1-exp(-lambda * (1-t))));
   CHIm_t = (delta_r)*(1-eta*(1-exp(-lambda * (1-t))));
   % Sigma_t= CHIm_t-CHIp_t;
   r_f_t  = (delta_r)*(1-eta);
   % r_f_t  = (1-eta)*Sigma_t;
end

gammad_t=lambda*sqrt(1/theta_t);
gammas_t=lambda*sqrt(theta_t);