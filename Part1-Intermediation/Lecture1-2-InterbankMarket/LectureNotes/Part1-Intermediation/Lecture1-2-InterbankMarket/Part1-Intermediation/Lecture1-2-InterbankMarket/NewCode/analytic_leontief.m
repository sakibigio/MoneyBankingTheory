function [CHIm_t,CHIp_t,r_f_t,gammad_t,gammas_t,theta_t]=analyticleontief(t,theta_0,delta_r,matchtech)
% alpha=matchtech.alpha;   % Technology for matching --- trade breakdown probability
lambda=matchtech.lambda; % Matching Efficienty
eta=matchtech.eta;       % Bargaining Power of Borrower
if theta_0 < 1
    theta_t = (1 + ((1-theta_0)/theta_0)*exp(lambda*t))^(-1);
    theta_1 = (1 + ((1-theta_0)/theta_0)*exp(lambda))^(-1);
    CHIp_t = (delta_r)*(theta_1/theta_t)^eta * ((theta_t-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    CHIm_t = (delta_r)*(theta_1/theta_t)^eta * ((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    r_f_t   = (theta_1/theta_t)^eta * (((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1))-eta*((1-theta_t)/(1-theta_1)))*(delta_r) ;
end
if theta_0 > 1
    theta_t = 1 + (theta_0-1)*exp(lambda*t);
    theta_1 = 1 + (theta_0-1)*exp(lambda);
    CHIp_t = (delta_r)*(theta_1/theta_t)^eta * ((theta_t-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    CHIm_t = (delta_r)*(theta_1/theta_t)^eta * ((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    r_f_t   = (theta_1/theta_t)^eta *(((theta_t^eta*theta_1^(1-eta)-1)/(theta_1-1))-eta*((theta_t-1)/(theta_1-1)))*(delta_r);
end
if theta_0 == 1
    theta_t = 1;
    CHIp_t = (delta_r) * (1-eta) *  (1-exp(-lambda*(1-t)));
    CHIm_t = (delta_r) * (1 - eta*(1-exp(-lambda*(1-t))));
    r_f_t   = (1-eta) * (delta_r) ;
end
gammad_t=lambda*min(1/theta_t,1);
gammas_t=lambda*min(theta_t,1)  ;