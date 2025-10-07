function [G]=build_matching(rho,lambda)
% rho is elascitiy of substituion
if rho==1
    G=@(a,b) lambda*a.^(1/2).*b.^(1/2);
elseif rho==0
    G=@(a,b) lambda*min(a,b)         ;
else
    G=@(a,b) lambda*(1/2)^(rho/(rho-1))*(a.^(1-1/rho)+b.^(1-1/rho)).^(rho/(rho-1));
end