function [CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t,CHIm_t,Sigma_t,t_vec]=interbanksolve_analyticcobbdouglas(delta_r,theta_0,N,matchtech)
%     PROJECT: Portfolio Theory with Settlement Frictions
%     AUTHORS: Javier Bianchi & Saki Bigio
%        CODE: Saki Bigio
%      UPDATE: 29/IV/2025
%      MODULE: Interbank Market Equilibrium Contingent to theta_0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleans and sets-up envirnoment
% SB: 04-25 changed it to make it look closer to model. 

% Initializing Variables
t_vec    = linspace(0, 1, N);
CHIm_t   = NaN(N,1);
CHIp_t   = NaN(N,1);
r_f_t    = NaN(N,1);
gammad_t = NaN(N,1);
gammas_t = NaN(N,1);
theta_t  = NaN(N,1);

% Evlution of Trading Probs
for tt=1:N
    [CHIm,CHIp,r_f,gammad,gammas,theta]=analytic_cobbdouglas(t_vec(tt),theta_0,delta_r,matchtech);
    theta_t(tt)=theta  ;
    gammas_t(tt)=gammas;
    gammad_t(tt)=gammad;
    CHIm_t(tt)=CHIm    ;
    CHIp_t(tt)=CHIp    ;
    r_f_t(tt)=r_f      ;
end
% Initial Yields
CHIm_0=CHIm_t(1);
CHIp_0=CHIp_t(1);

Sigma_t=CHIm_t-CHIp_t;
