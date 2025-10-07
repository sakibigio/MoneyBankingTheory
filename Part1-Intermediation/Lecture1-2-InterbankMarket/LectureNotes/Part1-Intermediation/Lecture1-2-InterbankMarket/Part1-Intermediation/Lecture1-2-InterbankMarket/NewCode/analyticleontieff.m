function [xi_d,xi_s,r_i,gammad,gammas,theta_t]=analyticleontief(delta_r,theta_o,N,matchtech)

if theta_0 < 1
    theta_t = (1 + ((1-theta_0)/theta_0)*exp(lambda*t))^(-1);
    theta_1 = (1 + ((1-theta_0)/theta_0)*exp(lambda))^(-1);
    CHIp_tl(i,j) = (iw-im)*(theta_1/theta_t)^eta * ((theta_t-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    CHIn_tl(i,j) = (iw-im)*(theta_1/theta_t)^eta * ((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    IF_tl(i,j)   = (theta_1/theta_t)^eta * (((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1))-eta*((1-theta_t)/(1-theta_1)))*(iw-im) + im;
end
if theta_0 > 1
    theta_t = 1 + (theta_0-1)*exp(lambda*t);
    theta_1 = 1 + (theta_0-1)*exp(lambda);
    CHIp_tl(i,j) = (iw-im)*(theta_1/theta_t)^eta * ((theta_t-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    CHIn_tl(i,j) = (iw-im)*(theta_1/theta_t)^eta * ((1-theta_t^eta*theta_1^(1-eta))/(1-theta_1));
    IF_tl(i,j)   = (theta_1/theta_t)^eta *(((theta_t^eta*theta_1^(1-eta)-1)/(theta_1-1))-eta*((theta_t-1)/(theta_1-1)))*(iw-im) + im;
end
if theta_0 == 1
    theta_t = 1;
    CHIp_tl(i,j) = (1-eta) * (iw-im) * (1-exp(-lambda*(1-t)));
    CHIn_tl(i,j) = (iw-im) * (1 - eta*(1-exp(-lambda*(1-t))));
    IF_tl(i,j)   = (1-eta) * (iw-im) + im;
end