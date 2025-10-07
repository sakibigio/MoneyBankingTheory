function [xi_d,xi_s,r_i,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_continuous(delta_r,theta_o,N,matchtech)
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
Delta=1/N;
% alpha=matchtech.alpha;   % Technology for matching --- trade breakdown probability
lambda=matchtech.lambda; % Matching Efficienty
eta=matchtech.eta       ;       % Bargaining Power of Borrower

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
    theta=theta_t(tt)                                                ;
    if theta>0&&theta<inf
        theta_t(tt+1)=theta*exp(-(G(1/theta,1)-G(theta,1)))          ;   
        theta=theta_t(tt)                                            ;
        gammas(tt)=G(1,theta)                                        ;
        gammad(tt)=G(1/theta,1)                                      ;
    elseif theta==0
        theta_t(tt+1)=0                                              ;   
        gammas(tt)=0                                                 ;
        gammad(tt)=NaN                                               ;
    elseif theta==inf
        theta_t(tt+1)=inf                                          ;   
        gammas(tt)=NaN                                             ;
        gammad(tt)=0                                               ;
    end
end

% Evolution of Valuations
if theta_t(end)>0
    for tt=N+1:-1:2
        r_i(tt-1)=(1-eta)*Vd(tt)+eta*Vs(tt);
        Vs(tt-1)=gammas(tt-1)*(r_i(tt-1))+(1-gammas(tt-1))*Vs(tt);
        Vd(tt-1)=gammad(tt-1)*(r_i(tt-1))+(1-gammad(tt-1))*Vd(tt);
    end
elseif theta_t(end)==0
    r_i(:)=0;
    Vs(:)=0;
    Vd(:)=0;
elseif theta_t(end)==inf
    r_i(:)=delta_r;
    Vs(:)=delta_r;
    Vd(:)=delta_r;
end

% Initial Yields
xi_d=Vd(1);
xi_s=Vs(1);

if isnan(xi_d)
    xi_d=delta_r;
end

if isnan(xi_s)
    xi_s=delta_r;
end
