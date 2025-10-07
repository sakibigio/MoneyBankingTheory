%% Interbank Market Solve
% (c) Saki Bigio
%
% Solve for the liquidity cost function in the interbank market model in
% "Portfolio Theory with Settlement Frictions".  
% Order of Plots
% by extreme case: Leontief, Cobb-Douglas
% plots in tau given theta, plots in theta, plots in lambda
clear; close all;
set(0,'defaulttextInterpreter','latex') 
format long;

% Select Folder
if strcmp(getenv('HOME'),'/Users/sakiclaudia')
    foldername='/Users/sakiclaudia/Dropbox/Settlements_portfolio/New Code/Figures/';
else
    foldername='/Users/sakibigio/Dropbox/Settlements_portfolio/New Code/Figures/';
end

% Color Ordering for Plots
newcolors = [0.00 0.00 1.0
             0.5 0.5 0.9
             0.25 0.75 0.9
             0.25 0.20 0.9
             0.7 0.7 0.7];

% Print Figure
plotit = 0;
printit= 1;

% Figure Scaling Properties
% Define base sizes
baseTick = 18;        % tick labels
baseAxes = 20;        % axes labels/title use multiplier
baseLegend = 20;      % legend text

% Set default fonts on all figures
set(groot, ...
    'DefaultLegendBox', 'off',...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',baseTick, ...
    'DefaultAxesLabelFontSizeMultiplier',baseAxes/baseTick, ...
    'DefaultAxesTitleFontSizeMultiplier',baseAxes/baseTick, ...
    'DefaultLegendFontName','Times','DefaultLegendFontSize',baseLegend, ...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultLineLineWidth', 4);

%% Parameters
% Interbank Market Tests
period=1    ; % Length of model period
rdw=0.0035  ; % Discount Window rate
rer=0.0025  ; % Rate on Excess Reserves
theta_o=0.8 ; 
N=10000     ; % # of rounds
BPS_scale=1e4*12;
T=1         ; % Trading Period Length
Delta=T/N   ;

% Plot Outcomes
round=1/N:1/N:1;

% To quarterly
rdw=rdw/period*BPS_scale    ; % Transform annual to period rate
rer=rer/period*BPS_scale    ;
delta_r=rdw-rer             ;

% Baseline Parameters
lambda_o= 1.2;
eta_o   = 0.5;

% Technology parameters
matchtech.rho = inf                        ; % inf for Leontieff matching
matchtech.lambda= lambda_o                      ; % Probability break down of interbank market --- not used in paper
% matchtech.poisson=0.5                    ; % Poisson Intensity of matches
% matchtech.alpha = matchtech.poisson*1*T/N; % correction term
matchtech.eta   = eta_o                      ; % Bargaining Power

% Final bargained rate
r_1=(1-matchtech.eta)*(rdw-rer);

% Vectors for Plots
N_theta=100;
ltheta_vec= linspace(-5,5,N_theta);
theta_vec = exp(ltheta_vec);
e_vec=theta_vec*0+1;

%% Comparison Various Matching Functions - Matching Rates
p_vec=[0 -1/4 -1 -2 -inf]; % p=(rho_vec-1)/(rho_vec);
legendCell = {'$p=0$ (Cobb-Douglas)', '$p=-1/4$', '$p=-1$ (Harmonic Mean)', '$p=-2$', '$p=\infty$ (Leontief)'};
rho_vec=1./(1-p_vec);
N_rho=length(p_vec);
ltheta_g_mat=zeros(N_rho,N_theta);
lambda=lambda_o;
for rr=1:N_rho
    rho=rho_vec(rr);
    G=build_matching(rho,lambda);
    ltheta_g_mat(rr,:)=-(G(theta_vec.^(-1),e_vec)'-G(theta_vec,e_vec)');
end
g_bound=-matchtech.lambda*(min(theta_vec.^(-1),e_vec)'-min(theta_vec,e_vec)');

figure('Name','Growth Comparison')
plot(ltheta_vec,ltheta_g_mat); hold on;
plot(ltheta_vec,g_bound,'k:','LineWidth',1.0);
hold on;
grid on; axis tight;
linestyleorder("mixedstyles")
colororder(newcolors);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{rate of change}');
legend(legendCell,'Location','northwest');

% xlabel('$\mathbf{ln(\theta_0)}$');
% ylabel('\textbf{rate of change}');
% legend(legendCell,'Location','northwest','Box','off');

if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_growthrate.pdf'],'BackgroundColor','none');
end

%% Theta Plots - Comparisons
matchtech.lambda=lambda_o;
N_theta          = 101;
theta0_bot =  log(0.25);
theta0_top =  log(1/0.25);
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
CHIp  = zeros(N_rho,length(ln_THETA));
CHIm  = zeros(N_rho,length(ln_THETA));
matchtech.lambda=lambda_o/N;
for jj = 1:N_theta
    theta_0 = exp(ln_THETA(jj));
    for rr=1:N_rho
        matchtech.rho=rho_vec(rr);
        [CHIm_0,CHIp_0,~,~,~,~]=interbanksolve_continuous(delta_r,theta_0,N,matchtech);
        CHIp(rr,jj) = CHIp_0;
        CHIm(rr,jj) = CHIm_0;
    end
end

figure('Name','Chi + (comparison)')
plot(ln_THETA,CHIp); hold on;
% plot(ltheta_vec,g_bound,'k:','LineWidth',1.0);
hold on;
grid on; axis tight;
linestyleorder("mixedstyles")
colororder(newcolors);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{BPS}');
% legend(legendCell,'Location','northwest','Box','off');
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_chip_comp.pdf'],'BackgroundColor','none');
end


figure('Name','Chi -(comparison)')
plot(ln_THETA,CHIm); hold on;
% plot(ltheta_vec,g_bound,'k:','LineWidth',1.0);
hold on;
grid on;  axis tight;
linestyleorder("mixedstyles")
colororder(newcolors);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{BPS}');
% legend(legendCell,'Location','northwest','Box','off');
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_chim_comp.pdf'],'BackgroundColor','none');
end


%% Comparison of Theta Evolution - Discrete Time
% This code solves the version of the model in discrete time:
theta_o=theta_o/5;
matchtech.lambda=lambda_o*Delta;
theta_t_mat=zeros(N_rho,N+1);
gammad_t_mat=zeros(N_rho,N+1);
gammas_t_mat=zeros(N_rho,N+1);
r_f_t_mat=zeros(N_rho,N);
chis_t_mat=zeros(N_rho,N+1);
chid_t_mat=zeros(N_rho,N+1);
for rr=1:N_rho
    rho=rho_vec(rr);
    matchtech.rho = rho;
    
    [~,~,r_t,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_discrete(delta_r,theta_o,N,matchtech);
    theta_t_mat(rr,:)=theta_t;
    gammad_t_mat(rr,:)=gammad;
    gammas_t_mat(rr,:)=gammas;
    r_f_t_mat(rr,:)=r_t;
    chis_t_mat(rr,:)=Vs;
    chid_t_mat(rr,:)=Vd;
end
matchtech.lambda=lambda_o;
[~,~,T_max]=probs_cobbdouglas(theta_o,matchtech);

figure('Name','Evolution Comparison (discrete time)')
plot([round 1],theta_t_mat); hold on;
hold on;
grid on; axis tight;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
ylabel('$\theta_{\tau}$')
linestyleorder("mixedstyles")
colororder(newcolors);
line([T_max T_max],[0 theta_o], 'LineStyle',':','Color','k');
text(T_max-0.02,theta_o/2.5,'End of Trade (Cobb Douglas case)', 'Interpreter', 'Latex', 'FontSize',12,'Rotation',90);
% legend(legendCell,'Location','northwest','Box','off');
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_trajectories.pdf'],'BackgroundColor','none');
end

if plotit==1
    figure('Name','Growth Comparison log (discrete time)')
    plot([round 1],log(theta_t_mat)); hold on;
    hold on;
    grid on; axis tight;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Matching Probs (deficit side) (discrete time)')
    plot([round 1],gammad_t_mat); hold on;
    hold on;
    grid on; axis tight;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$ ')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Matching Probs (surplus) (discrete time)')
    plot([round 1],gammas_t_mat); hold on;
    hold on;
    grid on; axis tight;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    
    figure('Name','Rates (discrete time)')
    plot(round,r_f_t_mat); hold on;
    plot(round,round*0+r_1,'LineWidth',1.0,'LineStyle',':'); hold on;
    ylim([0 rdw-rer]);
    hold on;
    grid on;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    
    figure('Name','Chi Surplus (discrete time)')
    plot([round 1],chis_t_mat); hold on;
    plot([round 1],[round 1]*0+r_1,'LineWidth',1.0,'LineStyle',':'); hold on;
    ylim([0 rdw-rer]);
    hold on;
    grid on;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Chi Deficit (discrete time)')
    plot([round 1],chid_t_mat); hold on;
    plot([round 1],[round 1]*0+r_1,'LineWidth',1.0,'LineStyle',':'); hold on;
    ylim([0 rdw-rer]);
    hold on;
    grid on;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
end

% Market Tightness
theta_dt_mat=theta_t_mat;

%% Comparison of Theta Evolution - Continuous Time Discretization
% This code solves the version of the model in discrete time:
% Update parameters
theta_t_mat=zeros(N_rho,N+1);
gammad_t_mat=zeros(N_rho,N+1);
gammas_t_mat=zeros(N_rho,N+1);
r_f_t_mat=zeros(N_rho,N);
chis_t_mat=zeros(N_rho,N+1);
chid_t_mat=zeros(N_rho,N+1);
for rr=1:N_rho
    rho=rho_vec(rr);
    matchtech.rho = rho;
    [~,~,r_t,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_continuous(delta_r,theta_o,N,matchtech);
    theta_t_mat(rr,:)=theta_t;
    gammad_t_mat(rr,:)=gammad;
    gammas_t_mat(rr,:)=gammas;
    r_f_t_mat(rr,:)=r_t;
    chis_t_mat(rr,:)=Vs;
    chid_t_mat(rr,:)=Vd;
end

if plotit==1
    figure('Name','Growth Comparison (continuous time)')
    plot([round 1],theta_t_mat); hold on;
    hold on;
    grid on; axis tight;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Growth Comparison log (continuous time)')
    plot([round 1],log(theta_t_mat)); hold on;
    hold on;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Matching Probs (continuous time)')
    plot([round 1],gammad_t_mat); hold on;
    hold on;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Matching Probs (continuous time)')
    plot([round 1],gammas_t_mat); hold on;
    hold on;
    grid on; axis tight;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Rates')
    plot(round,r_f_t_mat); hold on;
    plot(round,round*0+r_1,'LineWidth',1.0,'LineStyle',':'); hold on;
    ylim([0 rdw-rer]);
    hold on; grid on;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Chi Surplus  (continuous time)')
    plot([round 1],chis_t_mat); hold on;
    plot([round 1],[round 1]*0+r_1,'LineWidth',1.0,'LineStyle',':'); hold on;
    ylim([0 rdw-rer]);
    hold on; grid on;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Chi Deficit  (continuous time)')
    plot([round 1],chid_t_mat); hold on;
    plot([round 1],[round 1]*0+r_1,'LineWidth',1.0,'LineStyle',':'); hold on;
    ylim([0 rdw-rer]);
    hold on; grid on;
    xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
    % ylabel('\textbf{BPS}')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend(legendCell,'Location','northwest','Box','off');
    
    figure('Name','Discret vs. CT comparison')
    subplot(1,2,1);
    plot([round 1],theta_dt_mat); hold on;
    hold on;
    grid on; axis tight;
    linestyleorder("mixedstyles")
    colororder(newcolors); grid on;
    
    title('Discrete')
    subplot(1,2,2);
    plot([round 1],theta_t_mat); hold on;
    hold on;
    grid on; axis tight;
    linestyleorder("mixedstyles")
    colororder(newcolors); grid on;
    
    title('Continuous')
end

%% Comparison Between Discrete Time, Continuous Time and Analytical Formulas
% Leontief Case
[~,~,r_t,gammad,gammas,theta_t,Chip_t,Chim_t]=interbanksolve_discrete(delta_r,theta_o,N,matchtech);
theta_t_dt=theta_t;
r_t_dt=r_t;
gammad_t_dt=gammad;
gammas_t_dt=gammas;
Chip_t_dt=Chip_t;
Chim_t_dt=Chim_t;
[~,~,r_t,gammad,gammas,theta_t,Chip_t,Chim_t]=interbanksolve_continuous(delta_r,theta_o,N,matchtech);
theta_t_ct=theta_t;
r_t_ct=r_t;
gammad_t_ct=gammad;
gammas_t_ct=gammas;
Chip_t_ct=Chip_t;
Chim_t_ct=Chim_t;
matchtech.lambda=lambda_o;
[CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t,CHIm_t,Sigma_t,t_vec]=interbanksolve_analyticleontief(delta_r,theta_o,N,matchtech);

figure('Name','Comparisons (Leontieff Case)');
plot([round 1],theta_t_dt,'LineWidth',3); hold on;
plot([round 1],theta_t_ct,'LineWidth',3); hold on;
plot(t_vec,theta_t,'LineWidth',3); hold on;
legend('Discrete Time','Continuous Time','Analytic','Location','northwest','Box','off');
linestyleorder("mixedstyles")
colororder(newcolors);
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
ylabel('$\mathbf{\theta_{\tau}}$')
legend('Discrete Time','Continuous Time','Analytic','Location','northwest','Box','off');

figure('Name','Comparisons (Leontieff Case)');
plot([round],r_t_dt,'LineWidth',3); hold on;
plot([round],r_t_ct,'LineWidth',3); hold on;
plot(t_vec,r_f_t,'LineWidth',3); hold on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
ylabel('$\mathbf{\theta_{\tau}}$')
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
legend('Discrete Time','Continuous Time','Analytic','Location','northwest','Box','off');

% Cobb-Douglas Case
matchtech.rho = 1;
matchtech.lambda=lambda_o*Delta;
[~,~,r_t,gammad,gammas,theta_t,Chip_t,Chim_t]=interbanksolve_discrete(delta_r,theta_o,N,matchtech);
theta_t_dt=theta_t;
r_t_dt=r_t;
gammad_t_dt=gammad;
gammas_t_dt=gammas;
Chip_t_dt=Chip_t;
Chim_t_dt=Chim_t;
[~,~,r_t,gammad,gammas,theta_t,Chip_t,Chim_t]=interbanksolve_continuous(delta_r,theta_o,N,matchtech);
r_t_ct=r_t;
theta_t_ct=theta_t;
gammad_t_ct=gammad;
gammas_t_ct=gammas;
Chip_t_ct=Chip_t;
Chim_t_ct=Chim_t;
matchtech.lambda=lambda_o;
[CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t,CHIm_t,Sigma_t,t_vec]=interbanksolve_analyticcobbdouglas(delta_r,theta_o,N,matchtech);

figure('Name','Comparisons (Cobb-Douglas Case)');
plot([round 1],theta_t_dt,'LineWidth',3); hold on;
plot([round 1],theta_t_ct,'LineWidth',3,'LineStyle','--','Color','b'); hold on;
plot(t_vec,theta_t,'LineWidth',3,'LineStyle',':'); hold on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
ylabel('$\mathbf{\theta_{\tau}}$')
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
legend('Discrete Time','Continuous Time','Analytic','Location','northwest','Box','off');

figure('Name','Comparisons (Cobb-Douglas)');
plot(round,r_t_dt,'LineWidth',3); hold on;
plot(round,r_t_ct,'LineWidth',3,'LineStyle','--','Color','b'); hold on;
plot(t_vec,r_f_t,'LineWidth',3,'LineStyle',':'); hold on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
ylabel('$\mathbf{\r^{f}_{\tau}}$')
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
legend('Discrete Time','Continuous Time','Analytic','Location','northwest','Box','off');

figure('Name','Comparisons (Cobb-Douglas)');
plot([round 1],Chim_t_dt,'LineWidth',3); hold on;
plot([round 1],Chim_t_ct,'LineWidth',3,'LineStyle','--','Color','b'); hold on;
plot(t_vec,CHIm_t,'LineWidth',3,'LineStyle',':'); hold on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
ylabel('$\mathbf{\chi^{-}_{\tau}}$')
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
legend('Discrete Time','Continuous Time','Analytic','Location','northwest','Box','off');

figure('Name','Comparisons (Cobb-Douglas)');
plot([round 1],Chip_t_dt,'LineWidth',3); hold on;
plot([round 1],Chip_t_ct,'LineWidth',3,'LineStyle','--','Color','b'); hold on;
plot(t_vec,CHIp_t,'LineWidth',3,'LineStyle',':'); hold on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
ylabel('$\mathbf{\chi^{+}_{\tau}}$')
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
legend('Discrete Time','Continuous Time','Analytic','Location','northwest','Box','off');

%% Overall Checking Probabilities and Rates
matchtech.lambda=lambda_o;
theta_o=0.5;

% Check Leontief Case
[CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t,CHIm_t,Sigma_t,t_vec]=interbanksolve_analyticleontief(delta_r,theta_o,N,matchtech);
[Gammad_l,Gammas_l]=probs_leontief(theta_o,matchtech);
Gammad_l_test=1-exp(-sum(gammad_t)/N);
Gammas_l_test=1-exp(-sum(gammas_t)/N);
check_Gammad_l=Gammad_l_test-Gammad_l_test;
check_Gammas_l=Gammas_l_test-Gammas_l_test;
RF_0_l=CHIp_0/(Gammas_l);
RF_a_l=averagerate_leontief(theta_o,delta_r,matchtech);
check_rate=RF_0_l-RF_a_l;

% Check CobbDouglas Case
[CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t,CHIm_t,Sigma_t,t_vec]=interbanksolve_analyticcobbdouglas(delta_r,theta_o,N,matchtech);
[Gammad_cd,Gammas_cd]=probs_cobbdouglas(theta_o,matchtech);
Gammad_cd_test=1-exp(-sum(gammad_t)/N);
Gammas_cd_test=1-exp(-sum(gammas_t)/N);
check_Gammad_cd=Gammad_cd_test-Gammad_cd_test;
check_Gammas_cd=Gammas_cd_test-Gammas_cd_test;
RF_0_cd=CHIp_0/(Gammas_cd);
RF_a_cd=averagerate_cobbdouglas(theta_o,delta_r,matchtech);
check_rate=RF_0_cd-RF_a_cd;

%% Test of Dilation Property
matchtech.lambda=lambda_o;
theta_o=0.5;

% Check Leontief Case
cut=10/9;
[CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t,CHIm_t,Sigma_t,t_vec]=interbanksolve_analyticleontief(delta_r,theta_o,N,matchtech);
Chim_mid=Chim_t((N)/cut);
matchtech.lambda=lambda_o*(1-1/cut);
[CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t,CHIm_t,Sigma_t,t_vec]=interbanksolve_analyticleontief(delta_r,theta_t((N)/cut),N,matchtech);
test_dilation_l=CHIm_0-Chim_mid;

%% Cobb-Douglas - Time Plots
matchtech.lambda=lambda_o;
THETA_0 = [0.90, 1, 1/0.90]';
legendCell = {'$\theta=0.90$', '$\theta=1$', '$\theta=1/0.90$'};
CHIp_t  = zeros(N, length(THETA_0));
CHIm_t  = zeros(N, length(THETA_0));
RF_t    = zeros(N, length(THETA_0));
gammap_t  = zeros(N, length(THETA_0));
gammam_t  = zeros(N, length(THETA_0));
Sigma_t   = zeros(N, length(THETA_0));
for jj = 1:length(THETA_0)
    theta_0 = THETA_0(jj);
    [CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t_t,CHIm_t_t,Sigma,t_vec]=interbanksolve_analyticcobbdouglas(delta_r,theta_0,N,matchtech);
    CHIp_t(:,jj) = CHIp_t_t;
    CHIm_t(:,jj) = CHIm_t_t;
    RF_t(:,jj)   = r_f_t;
    gammap_t(:,jj) = gammas_t;
    gammam_t(:,jj) = gammad_t;
    Sigma_t(:,jj)   = Sigma;
end

% ****** Figures 1 ***** in the paper
figure('Name','Iterbank Rates')
plot(round, RF_t(:,1)); hold on;
plot(round, RF_t(:,2));
plot(round, RF_t(:,3));
ax = gca;
setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
% plot(TIME, IF_t(:,1), TIME, IF_t(:,2), TIME, IF_t(:,3))
% title('\textbf{INTERBANK INTEREST RATE}', 'Interpreter', 'Latex', 'FontSize',14)
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_InterbankRate_tau.pdf'],'BackgroundColor','none');
end

figure('Name','Chi +')
plot(round, CHIp_t(:,1)); hold on;
plot(round, CHIp_t(:,2));
plot(round, CHIp_t(:,3));
ax = gca;
setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_Chiplus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','Chi -')
plot(round, CHIm_t(:,1)); hold on;
plot(round, CHIm_t(:,2));
plot(round, CHIm_t(:,3));
ax = gca;
setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_Chiminus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','Surplus')
plot(round, Sigma_t(:,1)); hold on;
plot(round, Sigma_t(:,2));
plot(round, Sigma_t(:,3));
ax = gca;
setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_Surplus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','gamma +')
plot(round, gammap_t(:,1)); hold on;
plot(round, gammap_t(:,2));
plot(round, gammap_t(:,3));
ax = gca;
% setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
legend(legendCell,'Location','northwest')
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_gammaplus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','gamma -')
plot(round, gammam_t(:,1)); hold on;
plot(round, gammam_t(:,2));
plot(round, gammam_t(:,3));
ax = gca;
% setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')

if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_gammaminus_tau.pdf'],'BackgroundColor','none');
end

%% Cobb-Douglas - Theta Plots
matchtech.lambda=lambda_o;
N_theta          = 1001;
theta0_bot =  log(0.25);
theta0_top =  log(1/0.25);
matchtech.lambda=lambda_o;
% matchtech.lambda=2.44;
% N_theta          = 1001;
% theta0_bot =  log(0.5);
% theta0_top =  log(1/0.5);
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
CHIp  = zeros(length(ln_THETA),1);
CHIm  = zeros(length(ln_THETA),1);
RF    = zeros(length(ln_THETA),1);
T_conv = zeros(length(ln_THETA),1);
theta_bar= zeros(length(ln_THETA),1);
for jj = 1:length(ln_THETA)
    theta_0 = exp(ln_THETA(jj));
    [CHIm_0,CHIp_0,~,~,~,theta_1]=analytic_cobbdouglas(0,theta_0,delta_r,matchtech);
    [Gammam,Gammap,T_max]=probs_cobbdouglas(theta_0,matchtech);
    RF_0=CHIp_0/(Gammap);
    CHIp(jj) = CHIp_0;
    CHIm(jj) = CHIm_0;
    RF(jj)   = RF_0; % Fix. This is the unconditional one.
    T_conv(jj)=T_max;
    theta_bar(jj)=theta_1;
end

% Graphs the solutions
loc_min=find(T_conv==1,1,'first');
loc_max=find(T_conv==1,1,'last');

figure('Name','Theta Plots')
plot(ln_THETA, RF); hold on;
plot(ln_THETA, CHIm); 
plot(ln_THETA, CHIp);
line([0 0], [0 delta_r],'Color','k','LineWidth',1,'LineStyle',':');
axis([theta0_bot theta0_top -5 delta_r+5]);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{BPS}');
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
% legend('$\chi^-$', '$\chi^+$', '$\bar{r}^{f}$', 'Location', 'SouthEast','Box','Off','AutoUpdate','off', 'Interpreter', 'Latex');
legend('$r^{f}$', '$\chi^-$', '$\chi^+$', 'Interpreter','Latex', 'Location', 'SouthEast','Box','Off','AutoUpdate','off');
text(ln_THETA(80),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
text(ln_THETA(end-300),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
line([ln_THETA(loc_min) ln_THETA(loc_min)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
line([ln_THETA(loc_max) ln_THETA(loc_max)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
hold on; 
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_prices_theta.pdf'],'BackgroundColor','none');
end

figure('Name','Time of Market')
plot(ln_THETA, T_conv); hold on;
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('Final Trading Round');
% yticks([0]);
grid on;
% text(ln_THETA(80),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$');
line([ln_THETA(loc_min) ln_THETA(loc_min)], [min(T_conv) 1],'Color','k','LineStyle',':');
line([ln_THETA(loc_max) ln_THETA(loc_max)], [min(T_conv) 1],'Color','k','LineStyle',':');
hold on; 
linestyleorder("mixedstyles")
colororder(newcolors); grid on; axis tight;
% axis([theta0_bot theta0_top min(T_conv)-0.05 1+0.05]);
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_T_theta.pdf'],'BackgroundColor','none');
end

figure('Name','Theta_ bar')
index_aux=(loc_min:loc_max);
index_aux2=[(1:loc_min-1) (loc_max+1:N)];
plot(ln_THETA(index_aux),log(theta_bar(index_aux)),'LineWidth', 3);
grid on;
min_thetabar=min(log(theta_bar(index_aux))); 
max_thetabar=max(log(theta_bar(index_aux)));
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
text(ln_THETA(80),(1-matchtech.eta)*max_thetabar-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], 0*[max_thetabar max_thetabar],'Color','k','LineStyle','-','LineWidth',1);
line([ln_THETA(loc_min) ln_THETA(loc_min)], [min_thetabar max_thetabar],'Color','k','LineStyle',':','LineWidth',1);
line([ln_THETA(loc_max) ln_THETA(loc_max)], [min_thetabar max_thetabar],'Color','k','LineStyle',':','LineWidth',1);
hold on; 
text(ln_THETA(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
xlabel('$\mathbf{ln(\theta_0)}$')
ylabel('$\mathbf{ln(\theta_1)}$')
line([ln_THETA(1) ln_THETA(end)], [max_thetabar max_thetabar],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [min_thetabar min_thetabar],'Color','k','LineStyle','-','LineWidth',1);
line([0 0], [min_thetabar max_thetabar],'Color','k','LineWidth',1);
axis([ln_THETA(1) ln_THETA(end) min_thetabar max_thetabar]);
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_theta_bar.pdf'],'BackgroundColor','none');
end

% Cobb-Douglas - Dispersion Theta Plots
matchtech.lambda=lambda_o;
N_theta          = 1001;
N_mid      = (N + 1) / 2;
theta0_bot =  log(0.25);
theta0_top =  log(1/0.25);
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
% CHIp  = zeros(N,length(ln_THETA));
% CHIm  = zeros(N,length(ln_THETA));
RF_t    = zeros(N,length(ln_THETA));
Vol_t    = zeros(1,length(ln_THETA));
% T_conv= zeros(length(ln_THETA),1);
for jj = 1:length(ln_THETA)
    theta_0 = exp(ln_THETA(jj));
    [xi_d,xi_s,RF,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_analyticcobbdouglas(delta_r,theta_0,N,matchtech);
    RF_t(:,jj)=RF;
    [Gammam,Gammap,T_max]=probs_cobbdouglas(theta_0,matchtech);
    Vol_t(:,jj)=(1-Gammam)/Gammam;
end

figure('Name','Dispersion (Cobb Douglas)')
index_aux=(loc_min:loc_max);
index_aux2=[(1:loc_min-1) (loc_max+1:N)];
plot(ln_THETA(index_aux),max(RF_t(:,index_aux),[],1)-min(RF_t(:,index_aux),[],1),'LineWidth', 3);
grid on;
xlabel('$\mathbf{ln(\theta_0)}$')
ylabel('Dispersion Q (\textbf{BPS})'); 
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
text(ln_THETA(80),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
line([ln_THETA(loc_min) ln_THETA(loc_min)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
line([ln_THETA(loc_max) ln_THETA(loc_max)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
hold on; 
text(ln_THETA(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([0 0], [0 delta_r],'Color','k','LineWidth',1);
axis([ln_THETA(1) ln_THETA(end) -5 delta_r+5]);
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_Q_theta.pdf'],'BackgroundColor','none');
end

figure('Name','Volume Theta (Cobb Douglas)')
index_aux=(1:N);
index_aux2=[(1:loc_min-1) (loc_max+1:N)];
maxvol=max(Vol_t);
minvol=min(Vol_t);
plot(ln_THETA,Vol_t,'LineWidth', 3);
grid on;
xlabel('$\mathbf{ln(\theta_0)}$')
ylabel('Relative Volume I'); 
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
% text(ln_THETA(80),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
% line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
line([ln_THETA(loc_min) ln_THETA(loc_min)], [0 maxvol],'Color','k','LineStyle',':','LineWidth',1);
line([ln_THETA(loc_max) ln_THETA(loc_max)], [0 maxvol],'Color','k','LineStyle',':','LineWidth',1);
hold on; 
% text(ln_THETA(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], [maxvol maxvol],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([0 0], [0 delta_r],'Color','k','LineWidth',1);
axis([ln_THETA(1) ln_THETA(end) 0 maxvol]);
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_Vol_theta.pdf'],'BackgroundColor','none');
end

%% Cobb-Douglas - Eta-Theta Plots
N_theta          = 1001;
theta0_bot =  log(0.25);
theta0_top =  log(1/0.25);
matchtech.lambda=lambda_o;
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
eta_vec  = [0.25 0.5 0.75]; N_eta=3;
CHIp  = zeros(length(ln_THETA),N_eta);
CHIm  = zeros(length(ln_THETA),N_eta);
RF    = zeros(length(ln_THETA),N_eta);
T_conv = zeros(length(ln_THETA),N_eta);
theta_bar= zeros(length(ln_THETA),N_eta);
for jj = 1:length(ln_THETA)
        for ii=1:length(eta_vec)
            matchtech.eta=eta_vec(ii);
            theta_0 = exp(ln_THETA(jj));
            [CHIm_0,CHIp_0,~,~,~,theta_1]=analytic_cobbdouglas(0,theta_0,delta_r,matchtech);
            [Gammam,Gammap,T_max]=probs_cobbdouglas(theta_0,matchtech);
            RF_0=CHIp_0/(Gammap);
            CHIp(jj,ii) = CHIp_0;
            CHIm(jj,ii) = CHIm_0;
            RF(jj,ii)   = RF_0; % Fix. This is the unconditional one.
            T_conv(jj,ii)=T_max;
            theta_bar(jj,ii)=theta_1;
        end
end

% Graphs the solutions
loc_min=find(T_conv(:,1)==1,1,'first');
loc_max=find(T_conv(:,1)==1,1,'last');

figure('Name','Theta-Eta Plots')
% plot(ln_THETA, RF); hold on;
plot(ln_THETA, CHIm); hold on;
plot(ln_THETA, CHIp);
line([0 0], [0 delta_r],'Color','k','LineWidth',1,'LineStyle',':');
axis([theta0_bot theta0_top -5 delta_r+5]);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{BPS}');
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
% legend('$\chi^-$', '$\chi^+$', '$\bar{r}^{f}$', 'Location', 'SouthEast','Box','Off','AutoUpdate','off', 'Interpreter', 'Latex');
legend('$\chi^{-} (\eta=0.25)$','$\chi^{+} (\eta=0.25)$', '$\chi^{-} (\eta=0.5)$','$\chi^{+} (\eta=0.5)$', '$\chi^{-} (\eta=0.75)$','$\chi^{+} (\eta=0.75)$', 'Interpreter','Latex', 'Location', 'SouthEast','Box','Off','AutoUpdate','off');
text(ln_THETA(80),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
% line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
% text(ln_THETA(end-300),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
line([ln_THETA(loc_min) ln_THETA(loc_min)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
line([ln_THETA(loc_max) ln_THETA(loc_max)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
hold on; 
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_symmetry_theta.pdf'],'BackgroundColor','none');
end

%% Cobb-Douglas - Efficiency Plots 
matchtech.eta=0.5;
theta_0 = 0.75;
LAMBDA_vec = (0.5:0.01:2.8);
CHIp  = zeros(length(LAMBDA_vec),1);
CHIm  = zeros(length(LAMBDA_vec),1);
RF    = zeros(length(LAMBDA_vec),1);
T_conv= zeros(length(LAMBDA_vec),1);
Vol_t = T_conv;
for jj = 1:length(LAMBDA_vec)
    matchtech.lambda = LAMBDA_vec(jj);
    [CHIm_0,CHIp_0,~,~,~,~]=analytic_cobbdouglas(0,theta_0,delta_r,matchtech);
    [Gammam,Gammap,T_max]=probs_cobbdouglas(theta_0,matchtech);
    RF_0=CHIp_0/(Gammap);
    CHIp(jj) = CHIp_0;
    CHIm(jj) = CHIm_0;
    RF(jj)   = RF_0; % Fix. This is the unconditional one.
    T_conv(jj)=T_max;
    Vol_t(jj)=(1-Gammam)/Gammam;
end
% Graphs the solutions
loc_max=find(T_conv==1,1,'last');

figure('Name',"Lambda_Chi_theta125")
plot(LAMBDA_vec,CHIp(:)); hold on;
plot(LAMBDA_vec,CHIm(:)); axis tight;
plot(LAMBDA_vec,RF(:)); axis tight;
legend('$\chi^{-}$','$\chi^{+}$','$\bar{r}^f$','Box','off','Location','NorthEast','AutoUpdate','Off');
ylabel('\textbf{BPS}'); 
xlabel('$\mathbf{\lambda}$');
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
line([LAMBDA_vec(1) LAMBDA_vec(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
text(LAMBDA_vec(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([LAMBDA_vec(loc_max+1) LAMBDA_vec(loc_max+1)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
grid on; 
axis([LAMBDA_vec(1) LAMBDA_vec(end) -5 delta_r+5]);
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_prices_lambda.pdf'],'BackgroundColor','none');
end

if plotit==1
    figure('Name','Liq Premia (efficiency)')
    plot(LAMBDA_vec,CHIm(:)*0.8+CHIp(:)*0.2); hold on;
    % title('\textbf{Penalties}', 'Interpreter', 'Latex', 'FontSize',15)
    xlabel('$\mathbf{\lambda}$')
    ylabel('\textbf{Liquidity Premium (BPS)}')
    grid on;
    axis([LAMBDA_vec(1) LAMBDA_vec(end) -5 delta_r+5]);
    
    figure('Name','Funding Premia (efficiency)')
    plot(LAMBDA_vec,CHIm(:)*0.8-CHIp(:)*0.2); hold on;
    % title('\textbf{Penalties}', 'Interpreter', 'Latex', 'FontSize',15)
    xlabel('$\mathbf{\lambda}$')
    ylabel('\textbf{BPS}')
    grid on;
    axis([LAMBDA_vec(1) LAMBDA_vec(end) -50 delta_r+5]);
end

% Dispersion as function of Lambda
RF_t    = zeros(N,length(LAMBDA_vec));
for jj = 1:length(LAMBDA_vec)
    matchtech.lambda = LAMBDA_vec(jj);
    [xi_d,xi_s,RF,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_analyticcobbdouglas(delta_r,theta_0,N,matchtech);
    RF_t(:,jj)=RF;
end

figure('Name',"Lambda_Q")
index=(1:loc_max);
% plot(LAMBDA_vec(index),max(RF_t(:,index),[],1)-min(RF_t(:,index),[],1),'LineWidth',1); hold on;
plot(LAMBDA_vec(index),max(RF_t(:,index),[],1)-min(RF_t(:,index),[],1)); hold on;
plot(LAMBDA_vec(loc_max+1:end),max(RF_t(:,loc_max+1:end),[],1)-min(RF_t(:,loc_max+1:end),[],1)); hold on;
ylabel('Dispersion Q (\textbf{BPS})'); 
xlabel('$\mathbf{\lambda}$')
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
line([LAMBDA_vec(1) LAMBDA_vec(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
text(LAMBDA_vec(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([LAMBDA_vec(loc_max+1) LAMBDA_vec(loc_max+1)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
grid on;
axis([LAMBDA_vec(1) LAMBDA_vec(end) -5 delta_r+5]);
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_Q_lambda.pdf'],'BackgroundColor','none');
end

figure('Name','Volume Lambda (Cobb Douglas)')
max_vol=max(Vol_t);
% index=(1:loc_max);
% plot(LAMBDA_vec(index),max(RF_t(:,index),[],1)-min(RF_t(:,index),[],1),'LineWidth',1); hold on;
plot(LAMBDA_vec,Vol_t); hold on;
% plot(LAMBDA_vec(loc_max+1:end),Vol_t(loc_max+1:end)); hold on;
ylabel('Relative Volume I'); 
xlabel('$\mathbf{\lambda}$')
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
line([LAMBDA_vec(1) LAMBDA_vec(end)], [max_vol max_vol],'Color','k','LineWidth',1);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([LAMBDA_vec(loc_max+1) LAMBDA_vec(loc_max+1)], [0 max_vol],'Color','k','LineStyle',':','LineWidth',1);
grid on;
axis([LAMBDA_vec(1) LAMBDA_vec(end) 0 max_vol]); 
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_cd_Vol_lambda.pdf'],'BackgroundColor','none');
end

%% Leontief - Time Plots
matchtech.lambda=lambda_o;
THETA_0 = [0.90, 1, 1/0.90]';
legendCell = {'$\theta=0.90$', '$\theta=1$', '$\theta=1/0.90$'};
CHIp_t  = zeros(N, length(THETA_0));
CHIm_t  = zeros(N, length(THETA_0));
RF_t    = zeros(N, length(THETA_0));
gammap_t  = zeros(N, length(THETA_0));
gammam_t  = zeros(N, length(THETA_0));
Sigma_t   = zeros(N, length(THETA_0));
for jj = 1:length(THETA_0)
    theta_0 = THETA_0(jj);
    [CHIm_0,CHIp_0,r_f_t,gammad_t,gammas_t,theta_t,CHIp_t_t,CHIm_t_t,Sigma,t_vec]=interbanksolve_analyticleontief(delta_r,theta_0,N,matchtech);
    CHIp_t(:,jj) = CHIp_t_t;
    CHIm_t(:,jj) = CHIm_t_t;
    RF_t(:,jj)   = r_f_t;
    gammap_t(:,jj) = gammas_t;
    gammam_t(:,jj) = gammad_t;
    Sigma_t(:,jj)   = Sigma;
end

figure('Name','Iterbank Rates')
plot(round, RF_t(:,1)); hold on;
plot(round, RF_t(:,2));
plot(round, RF_t(:,3));
ax = gca;
setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
grid on;
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_InterbankRate_tau.pdf'],'BackgroundColor','none');
end

figure('Name','Chi +')
plot(round, CHIp_t(:,1)); hold on;
plot(round, CHIp_t(:,2));
plot(round, CHIp_t(:,3));
ax = gca;
setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_Chiplus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','Chi -')
plot(round, CHIm_t(:,1)); hold on;
plot(round, CHIm_t(:,2));
plot(round, CHIm_t(:,3));
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_Chiminus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','Surplus')
plot(round, Sigma_t(:,1)); hold on;
plot(round, Sigma_t(:,2));
plot(round, Sigma_t(:,3));
ax = gca;
setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on; 
ylabel('\textbf{BPS}')
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_Surplus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','gamma +')
plot(round, gammap_t(:,1)); hold on;
plot(round, gammap_t(:,2));
plot(round, gammap_t(:,3));
linestyleorder("mixedstyles")
grid on;
colororder(newcolors); grid on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
legend(legendCell,'Location','northwest');
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_gammaplus_tau.pdf'],'BackgroundColor','none');
end

figure('Name','gamma -')
plot(round, gammam_t(:,1)); hold on;
plot(round, gammam_t(:,2));
plot(round, gammam_t(:,3));
ax = gca;
% setNumYTicks(ax, 5);
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
xlabel('\textbf{Trading Round} $\mathbf{\tau}$')
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_gammaminus_tau.pdf'],'BackgroundColor','none');
end


%% Leontief - Theta Plots
matchtech.lambda=lambda_o;
N_theta          = 1001;
N_mid      = (N + 1) / 2;
theta0_bot =  log(0.05);
theta0_top =  log(1/0.05);
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
CHIp  = zeros(length(ln_THETA),1);
CHIm  = zeros(length(ln_THETA),1);
RF    = zeros(length(ln_THETA),1);
T_conv= zeros(length(ln_THETA),1);
for jj = 1:length(ln_THETA)
    theta_0 = exp(ln_THETA(jj));
    [CHIm_0,CHIp_0,~,~,~,~]=analytic_leontief(0,theta_0,delta_r,matchtech);
    [Gammam,Gammap]=probs_leontief(theta_0,matchtech);
    RF_0=CHIp_0/(Gammap);
    CHIp(jj) = CHIp_0;
    CHIm(jj) = CHIm_0;
    RF(jj)   = RF_0; % Fix. This is the unconditional one.
    T_conv(jj)=T_max;
end

% Graphs the solutions
figure('Name','Theta Plots')
plot(ln_THETA, RF); hold on;
plot(ln_THETA, CHIm); 
plot(ln_THETA, CHIp);
line([0 0], [0 delta_r],'Color','k','LineWidth',1,'LineStyle',':');
axis([theta0_bot theta0_top -5 delta_r+5]);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{BPS}');
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
legend('$r^{f}$', '$\chi^-$', '$\chi^+$', 'Interpreter','Latex', 'Location', 'SouthEast','Box','Off','AutoUpdate','off');
text(ln_THETA(80),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
text(ln_THETA(end-300),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
% line([ln_THETA(loc_min) ln_THETA(loc_min)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
% line([ln_THETA(loc_max) ln_THETA(loc_max)], [0 delta_r],'Color','k','LineStyle',':','LineWidth',1);
hold on; 
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_prices_theta.pdf'],'BackgroundColor','none');
end

% Dispersion plot
matchtech.lambda=lambda_o;
N_theta          = 1001;
N_mid      = (N + 1) / 2;
theta0_bot =  log(0.25);
theta0_top =  log(1/0.25);
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
RF_t    = zeros(N,length(ln_THETA));
Vol_t   = zeros(1,length(ln_THETA));
for jj = 1:length(ln_THETA)
    theta_0 = exp(ln_THETA(jj));
    [xi_d,xi_s,RF,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_analyticleontief(delta_r,theta_0,N,matchtech);
    RF_t(:,jj)=RF;
    [Gammam,Gammap]=probs_leontief(theta_0,matchtech);
    Vol_t(jj)=(1-Gammam)/Gammam;
end

figure('Name','Dispersion (Leontief)')
plot(ln_THETA,max(RF_t,[],1)-min(RF_t,[],1),'LineWidth', 3);
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
text(ln_THETA(80),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
hold on; 
text(ln_THETA(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
xlabel('$\mathbf{ln(\theta_0)}$')
ylabel('Dispersion Q (\textbf{BPS})'); 
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([0 0], [0 delta_r],'Color','k','LineWidth',1);
axis([ln_THETA(1) ln_THETA(end) -5 delta_r+5]);
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_Q_theta.pdf'],'BackgroundColor','none');
end

figure('Name','Volume Theta (Leontief)')
maxvol=max(Vol_t);
minvol=min(Vol_t);
plot(ln_THETA,Vol_t,'LineWidth', 3);
grid on;
xlabel('$\mathbf{ln(\theta_0)}$')
ylabel('Relative Volume I'); 
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
hold on; 
line([ln_THETA(1) ln_THETA(end)], [maxvol maxvol],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([0 0], [0 delta_r],'Color','k','LineWidth',1);
axis([ln_THETA(1) ln_THETA(end) 0 maxvol]);
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_Vol_theta.pdf'],'BackgroundColor','none');
end

%% Leontief - Eta-Theta Plots
matchtech.lambda=lambda_o/1.5;
N_theta          = 1001;
theta0_bot =  log(0.05);
theta0_top =  log(1/0.05);
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
eta_vec  = [0.25 0.5 0.75]; N_eta=3;
CHIp  = zeros(length(ln_THETA),N_eta);
CHIm  = zeros(length(ln_THETA),N_eta);
RF    = zeros(length(ln_THETA),N_eta);
T_conv = zeros(length(ln_THETA),N_eta);
theta_bar= zeros(length(ln_THETA),N_eta);
for jj = 1:length(ln_THETA)
        for ii=1:length(eta_vec)
            matchtech.eta=eta_vec(ii);
            theta_0 = exp(ln_THETA(jj));
            [CHIm_0,CHIp_0,~,~,~,~]=analytic_leontief(0,theta_0,delta_r,matchtech);
            CHIp(jj,ii) = CHIp_0;
            CHIm(jj,ii) = CHIm_0;            
        end
end

figure('Name','Theta-Eta Plots')
% plot(ln_THETA, RF); hold on;
plot(ln_THETA, CHIm); hold on;
plot(ln_THETA, CHIp);
line([0 0], [0 delta_r],'Color','k','LineWidth',1,'LineStyle',':');
axis([theta0_bot theta0_top -5 delta_r+5]);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{BPS}');
grid on;
linestyleorder("mixedstyles")
legend('$\chi^{-} (\eta=0.25)$','$\chi^{+} (\eta=0.25)$', '$\chi^{-} (\eta=0.5)$','$\chi^{+} (\eta=0.5)$', '$\chi^{-} (\eta=0.75)$','$\chi^{+} (\eta=0.75)$', 'Interpreter','Latex', 'Location', 'East','Box','Off','AutoUpdate','off');
colororder(newcolors); grid on;
% text(ln_THETA(80),(1-matchtech.eta)*delta_r-6,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
% line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-','LineWidth',1);
hold on; 
text(ln_THETA(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
xlabel('$\mathbf{ln(\theta_0)}$')
ylabel('\textbf{BPS}')
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([ln_THETA(1) ln_THETA(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
line([0 0], [0 delta_r],'Color','k','LineWidth',1);
axis([ln_THETA(1) ln_THETA(end) -5 delta_r+5]);
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_symmetry_theta.pdf'],'BackgroundColor','none');
end

%% Leontief - Efficiency Plots 
theta_0 = 0.75;
matchtech.lambda=lambda_o;
matchtech.eta=0.5;
LAMBDA_vec = (0.5:0.1:3);
CHIp  = zeros(length(LAMBDA_vec),1);
CHIm  = zeros(length(LAMBDA_vec),1);
RF    = zeros(length(LAMBDA_vec),1);
Vol_t = zeros(length(LAMBDA_vec),1);
for jj = 1:length(LAMBDA_vec)
    matchtech.lambda = LAMBDA_vec(jj);
    [CHIm_0,CHIp_0,~,~,~,~]=analytic_leontief(0,theta_0,delta_r,matchtech);
    [Gammam,Gammap]=probs_leontief(theta_0,matchtech);
    RF_0=CHIp_0/(Gammap);
    CHIp(jj) = CHIp_0;
    CHIm(jj) = CHIm_0;
    RF(jj)   = RF_0; % Fix. This is the unconditional one.
    Vol_t(jj)= (1-Gammap)/Gammap;
end

figure('Name',"Prices Lambda (Loentief)")
plot(LAMBDA_vec,CHIp(:)); hold on;
plot(LAMBDA_vec,CHIm(:)); axis tight;
plot(LAMBDA_vec,RF(:)); axis tight;
legend('$\chi^{-}$','$\chi^{+}$','$\bar{r}^f$','Box','off','Location','SouthEast','AutoUpdate','Off');
ylabel('\textbf{BPS}'); 
xlabel('$\mathbf{\lambda}$');
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
line([LAMBDA_vec(1) LAMBDA_vec(end)], [delta_r delta_r],'Color','k','LineWidth',1);
text(LAMBDA_vec(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
grid on; 
axis([LAMBDA_vec(1) LAMBDA_vec(end) -5 delta_r+5]);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_prices_lambda.pdf'],'BackgroundColor','none');
end

if plotit==1
    figure('Name','Liq Premia (efficiency)')
    plot(LAMBDA_vec,CHIm(:)*0.8+CHIp(:)*0.2); hold on;
    % title('\textbf{Penalties}', 'Interpreter', 'Latex', 'FontSize',15)
    xlabel('$\mathbf{\lambda}$')
    ylabel('\textbf{Liquidity Premium (BPS)}')
    grid on;
    axis([LAMBDA_vec(1) LAMBDA_vec(end) -5 delta_r+5]);
    
    figure('Name','Funding Premia (efficiency)')
    plot(LAMBDA_vec,CHIm(:)*0.8-CHIp(:)*0.2); hold on;
    % title('\textbf{Penalties}', 'Interpreter', 'Latex', 'FontSize',15)
    xlabel('$\mathbf{\lambda}$')
    ylabel('\textbf{BPS}')
    grid on;
    axis([LAMBDA_vec(1) LAMBDA_vec(end) -50 delta_r+5]);
end

% Dispersion as function of Lambda
theta_0 = 0.75;
RF_t    = zeros(N,length(LAMBDA_vec));
for jj = 1:length(LAMBDA_vec)
    matchtech.lambda = LAMBDA_vec(jj);
    [xi_d,xi_s,RF,gammad,gammas,theta_t,Vs,Vd]=interbanksolve_analyticleontief(delta_r,theta_0,N,matchtech);
    RF_t(:,jj)=RF;
end

figure('Name','Lambda_Q (Leontief)')
plot(LAMBDA_vec,max(RF_t,[],1)-min(RF_t,[],1)); hold on;
% title('\textbf{Max-Min Interbank Rates}', 'Interpreter', 'Latex', 'FontSize',15)
ylabel('\textbf{BPS}')
xlabel('$\mathbf{\lambda}$')
grid on;
linestyleorder("mixedstyles")
ylabel('Dispersion Q (\textbf{BPS})'); 
xlabel('$\mathbf{\lambda}$');
grid on;
linestyleorder("mixedstyles")
axis([LAMBDA_vec(1) LAMBDA_vec(end) -5 delta_r+5]);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
colororder(newcolors); grid on;
line([LAMBDA_vec(1) LAMBDA_vec(end)], [delta_r delta_r],'Color','k','LineWidth',1);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
text(LAMBDA_vec(5),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
axis tight;
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_Q_lambda.pdf'],'BackgroundColor','none');
end

figure('Name','Volume Lambda (Leontief)')
max_vol=max(Vol_t);
% index=(1:loc_max);
% plot(LAMBDA_vec(index),max(RF_t(:,index),[],1)-min(RF_t(:,index),[],1),'LineWidth',1); hold on;
plot(LAMBDA_vec,Vol_t); hold on;
% plot(LAMBDA_vec(loc_max+1:end),Vol_t(loc_max+1:end)); hold on;
ylabel('Relative Volume I'); 
xlabel('$\mathbf{\lambda}$')
grid on;
linestyleorder("mixedstyles")
colororder(newcolors); grid on;
line([LAMBDA_vec(1) LAMBDA_vec(end)], [max_vol max_vol],'Color','k','LineWidth',1);
line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 0],'Color','k','LineStyle','-','LineWidth',1);
% line([LAMBDA_vec(1) LAMBDA_vec(end)], [0 max_vol],'Color','k','LineStyle',':','LineWidth',1);
grid on;
axis([LAMBDA_vec(1) LAMBDA_vec(end) 0 max_vol]);
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[foldername 'F_l_Vol_lambda.pdf'],'BackgroundColor','none');
end

%% Walrasian Limit
matchtech.lambda=lambda_o*100;
N_theta          = 1001;
theta0_bot =  log(0.05);
theta0_top =  log(1/0.05);
ln_THETA = linspace(theta0_bot, theta0_top, N_theta);
CHIp  = zeros(N_theta,1);
CHIm  = zeros(N_theta,1);
RF    = zeros(N_theta,1);
T_conv= zeros(N_theta,1);
for jj = 1:N_theta
    theta_0 = exp(ln_THETA(jj));
    [CHIm_0,CHIp_0,~,~,~,~]=analytic_leontief(0,theta_0,delta_r,matchtech);
    [Gammam,Gammap]=probs_leontief(theta_0,matchtech);
    RF_0=CHIp_0/(Gammap);
    CHIp(jj) = CHIp_0;
    CHIm(jj) = CHIm_0;
    RF(jj)   = RF_0; % Fix. This is the unconditional one.
    T_conv(jj)=T_max;
end

figure('Name','Theta Plots')
N_mid=((N_theta+1)/2);
index1=(1:N_mid-1);
index2=(N_mid+1:N_theta);
plot(ln_THETA(index1), RF(index1)','Color',newcolors(1,:)); hold on;
plot(ln_THETA(index2), RF(index2)','Color',newcolors(1,:)); hold on;
scatter(ln_THETA(N_mid), (1-matchtech.eta)*delta_r,40,'MarkerFaceColor',newcolors(1,:),'MarkerEdgeColor',newcolors(1,:));
scatter(ln_THETA(N_mid), 0*delta_r,40,'MarkerFaceColor','w','MarkerEdgeColor',newcolors(1,:));
scatter(ln_THETA(N_mid), delta_r,40,'MarkerFaceColor','w','MarkerEdgeColor',newcolors(1,:));
line([0 0], [0 delta_r],'Color','k', 'LineWidth',2,'LineStyle',':');
axis([theta0_bot theta0_top -5 delta_r+5]);
xlabel('$\mathbf{ln(\theta_0)}$');
ylabel('\textbf{BPS}');
grid on;

legend('$\chi^-=\chi^+=r^{f}$', 'Location', 'best','Box','Off','AutoUpdate','off','Interpreter','Latex');
text(ln_THETA(80),delta_r-6,'$\mathbf{r^{w}-r^{m}}$','FontSize',16);
line([ln_THETA(1) ln_THETA(end)], [delta_r delta_r],'Color','k', 'LineWidth',2);
line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[delta_r delta_r],'Color','k','LineStyle','-', 'LineWidth',2);
%line([ln_THETA(1) ln_THETA(end)], (1-exp(-matchtech.eta*lambda_o))*[delta_r delta_r],'Color','k','LineStyle',':');
%text(ln_THETA(600),(1-exp(-matchtech.eta*lambda_o))*delta_r+3,'$(1-\exp(-\eta\bar{\lambda}))*(\mathbf{r^{w}-r^{m}})$');
%line([ln_THETA(1) ln_THETA(end)], (exp(-(1-matchtech.eta)*lambda_o))*[delta_r delta_r],'Color','k','LineStyle',':');
%text(ln_THETA(600),exp(-(1-matchtech.eta)*lambda_o)*delta_r+3,'$\exp(-(1-\eta)\bar{\lambda})*(\mathbf{r^{w}-r^{m}})$');
line([ln_THETA(1) ln_THETA(end)], (1-matchtech.eta)*[0 0],'Color','k','LineStyle','-', 'LineWidth',2);
text(ln_THETA(80),(1-matchtech.eta)*delta_r+3,'$(1-\eta)*(\mathbf{r^{w}-r^{m}})$','FontSize',16);
hold on; 
if printit==1
    orient landscape;
    % saveas(gcf, 'Dist_example', 'pdf')
    ax = gca;
    exportgraphics(ax,[foldername 'F_Walrasian_prices_theta.pdf'],'BackgroundColor','none');
end


%% Functions for plots
function setNumYTicks(ax, N)
    limits = ax.YLim;
    ax.YTick = linspace(limits(1), limits(2), N);
end