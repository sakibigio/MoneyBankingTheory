function [xi_d,xi_s,r_i,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_discrete(delta_r,theta_o,N,matchtech)
%     PROJECT: Portfolio Theory with Settlement Frictions
%     AUTHORS: Javier Bianchi & Saki Bigio
%        CODE: Saki Bigio
%      UPDATE: 29/IV/2025
%      MODULE: Interbank Market Equilibrium Contingent to theta_0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleans and sets-up envirnoment
% SB: 04-25 changed it to make it look closer to model. 

% Bring Parameters
rho=matchtech.rho;       % Technology for matching --- CES p=1-1/rho-> rho=0 Leontieff, rho=1/2, cobb douglas=1. 
% alpha=matchtech.alpha;   % Technology for matching --- trade breakdown probability
lambda=matchtech.lambda; % Matching Efficienty
eta=matchtech.eta;       % Bargaining Power of Borrower

% Construct Matching Function - with Efficiency embedded
G=build_matching(rho,lambda);

% Initializing Variables
theta_t = theta_o*ones(N,1);
r_i     = ones(N,1)        ;
gammad  = zeros(N+1,1)     ;
gammas  = zeros(N+1,1)     ;

% New Variable
Vd    = delta_r*ones(N+1,1)  ;
Vs    = zeros(N+1,1)  ;

% Evlution of Trading Probs
for tt=1:N
    theta=theta_t(tt)                                            ;
    theta_t(tt+1)=theta*max(1-G(1/theta,1),0)/max(1-G(1,theta),0);   
    theta=theta_t(tt)                                            ;
    gammas(tt)=min(G(1,theta),1)                                        ;
    gammad(tt)=min(G(1/theta,1),1)                                      ;
end

% Evolution of Valuations
for tt=N+1:-1:2
    r_i(tt-1)=(1-eta)*Vd(tt)+eta*Vs(tt);
    Vs(tt-1)=gammas(tt-1)*(r_i(tt-1))+(1-gammas(tt-1))*Vs(tt);
    Vd(tt-1)=gammad(tt-1)*(r_i(tt-1))+(1-gammad(tt-1))*Vd(tt);
end

% Initial Yields
xi_d=Vd(1);
xi_s=Vs(1);


