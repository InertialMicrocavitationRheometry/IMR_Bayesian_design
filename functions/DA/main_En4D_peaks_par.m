function [x,E2,P_post,params]= main_En4D_peaks_par(Model_params,t,yth,R0_all,Req_all,tspan_all,peak_time_idx)

% parallel En4D_Var Ensemble Kalman filter for IMR solve

% based on Sakov 2011 paper and Spratt et al. 2021 paper
% Original author: Jean-Sebastien Spratt (jspratt@caltech.edu)
% Edited by: Tianyi Chu (tchu72@gatech.edu)


%% 

q               =    48;  % Number of ensemble members
exp_i           =    1;


%% Data assimilation parameters

method          =   'En4D'; % data assimilation method ('En4D','EnKS',('EnKF'))

% Initial parameter guesses (all must be specified even if not used, in
% order to run):
model           =    Model_params{1};
P_prior         =    Model_params{3};
mu_theta        =    P_prior.mu; 
sigma_theta     =    P_prior.sigma;
theta_params    =    mvnrnd(mu_theta,sigma_theta,q);
G_prior         =    theta_params(:,1) ;
mu_prior        =    theta_params(:,2)  ;
alpha_prior     =    theta_params(:,3);
G_guess         =    mean(G_prior);
mu_guess        =    mean(mu_prior);
alpha_guess     =    mean(alpha_prior);
Input_prior     =    true;
G1_guess        =    1e9;
lambda_nu_guess =    0.1;
init_scheme     =    2; % leave as 2, initializes ensemble with truth + noise

% The following are ending criteria for iterative optimization:
epsilon         =    1e-10; % threshold difference in norm of state vector between steps
max_iter        =    5; % max # of iterations until end optimization (5 to 10 is usually good)

% Note: Beta coefficients are fixed as equal here. They can be modified in
% the main_En4D_peaks and main_mda_peaks files if needed, if more weight
% needs to be attributed to earlier or later points in assimilation

%% Modeling parameters

NT = 240; % Amount of nodes inside the bubble
NTM = 240; % Amount of nodes outside the bubble
Pext_type = 'IC'; % Type of external forcing
ST = 0.056; % (N/m) Liquid Surface Tension
Tgrad = 1; % Thermal effects inside bubble
Tmgrad = 1; % Thermal effects outside bubble
Cgrad = 1; % Vapor diffusion effects
comp = 1; % Activates the effect of compressibility (0=Rayleigh-Plesset, 1=Keller-Miksis)
disp_timesteps = 1; % displays timesteps in En4D run (for debugging)

% following should not be changed (untested):
disptime = 0; % 1 = display simulation time
Dim = 0; % 1 = output variables in dimensional form



%% Covariance Inflation parameters
% The following is only used for the IEnKS. Default values provided below
% See Spratt et al. (2020) section 3.3 for details on theta, lambda

CI_scheme = 2; % 1 is scalar alpha CI, 2 is RTPS (BEST)
CI_theta = 0.7; % multiplicative covariance parameter (0.5 < theta < 0.95)
CI_add = 0; % Set to 1 for additive covariance (else 0)

beta = 1.02; % additive covariance parameter (lambda in paper) (1.005 < beta < 1.05)

alpha = 0.005; % random noise in forecast step (only for EnKF)


%% Spread of parameters in the ensemble
%
%Rspread = 0.02;
Rspread = 0.00;
Uspread = 0.0;
Pspread = 0.0;
Sspread = 0.0;
tauspread = 0.0;
Cspread = 0.000;
Tmspread = 0.0000;
Brspread = 0.00;
fohspread = 0.00;
Despread = 0; % set to 0 if not used in model
lambda_nuspread = 0; % set to 0 if not used in model



%% Do not modify
visco_params = struct('G',G_guess,'G1',G1_guess,'mu',mu_guess, ...
    'alpha',alpha_guess,'lambda_nu',lambda_nu_guess);
est_params = [];


%%

%rng(99)
tic

%% Initialize and import data

clear aux1 aux2 dy2

n =  length(t)-1;

% Pre-allocating memory:

x          =    zeros(2*NT+NTM+11,q,n+1);
x_est      =    zeros(2*NT+NTM+11,n+1);
E_est      =    zeros(2*NT+NTM+11,q,max_iter);
E2         =    zeros(q,max(peak_time_idx)-1);

tspan      =    tspan_all(exp_i);
R0         =    R0_all(exp_i);
Req        =    Req_all(exp_i);

create_ensemble;

tspan_star =    max(tspan);
x(:,:,1)   =    xi;
tau_del    =    cell(q,1);
idx        =    1;


%% Iterate

vars = {NT Pext_type Pext_Amp_Freq disptime Tgrad Tmgrad ...
    Cgrad comp t0 neoHook nhzen sls linkv k chi fom foh We Br A_star ...
    B_star Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star ...
    De deltaY yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm tspan_star NTM rho R0 fung fung2 fungexp fungnlvis};

opts.POSDEF = true;
opts.SYM = true;


%delta = 1.005; %inflation param
delta = beta;
%epsilon = 0.001;


% Special case: deleting points around first and second collapse
%{
num_del_pts = 1; % num points to delete around collapse (can be 0 for 
                 % only collapse point deleted)
if num_peaks > 2
    collapse_1_idx = peak_indices(2)-1;
else
    collapse_1_idx = peak_indices(1)-1;
end
Beta(collapse_1_idx-num_del_pts:collapse_1_idx+num_del_pts) = 0;
Beta(l-num_del_pts:l) = 0;
%}


l                   =   max(peak_time_idx(:)) - 1; % size of data assimilation window
Beta                =   (1/l).*ones(1,l); % MDA coefs (equal here)
timestep_time       =   zeros(1,n-l+1);
time_index          =   1;
j_max               =   3;                % Numbers of DA cycles
params              =   [];

for j = 1:j_max

    %  Restrict mu and G to be positive 

    if min(x(2*NT+NTM+7,:,j))<=0
        Idx = (x(2*NT+NTM+7,:,j))<=0;
        x(2*NT+NTM+7,Idx,j) = x(2*NT+NTM+7,Idx,j)-min(x(2*NT+NTM+7,:,j))+1e-6;
    end
    if min(x(2*NT+NTM+8,:,j))<=0
        Idx = (x(2*NT+NTM+8,:,j))<=0;
        x(2*NT+NTM+8,Idx,j) = x(2*NT+NTM+8,Idx,j)-min(x(2*NT+NTM+8,:,j))+1e-6;
    end

    %j % print index for debugging if it crashes
    %timestep_time(j) = toc;
    x10 =  mean(squeeze(x(:,:,j)),2);
    A10 =  squeeze(x(:,:,j)) - x10*ones(1,q);
    x1  =  x10;
    TT  =  eye(q);
    dx1 =  1; %initialize greater than epsilon
    jj  =  1;
    A1  =  A10*TT;
    E1  =  x1*ones(1,q) + A1;

    for j1 = 1:size(A1,1)
        if std(A1(j1,:)) ==0
            std_prior(j1) =0;
        else
            std_prior(j1) = sqrt(robustcov(A1(j1,:)));
        end
    end


    while norm(dx1) > epsilon & jj <= max_iter


        %%%%%%%%%%%%%%-----saving parameters at each iteration-----%%%%%%%%%%%%%%%%%%

        G_all            =  E1(2*NT+NTM+7,:)*P_inf;
        mu_all           =  E1(2*NT+NTM+8,:)/Uc*(P_inf*R0);
        alpha_all        =  (E1(2*NT+NTM+10,:));
        params(:,1,idx)  =  G_all;
        params(:,2,idx)  =  mu_all;
        params(:,3,idx)  =  alpha_all;
        idx              =  idx+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        TTinv = linsolve(TT,eye(size(TT)));
        t1    = t(exp_i,time_index);
        t2    = t(exp_i,time_index+l);

        parfor memb = 1:q

            [t_memb{memb}, EE{memb},~] =  f_new(t1,t2,E1(:,memb),vars,tau_del{memb});
            t_sim                      =  t_memb{memb};
            y_sim                      =  EE{memb}(:,1);
            U_sim                      =  EE{memb}(:,2);
            y2(1,memb,:)               =  interp1(t_sim,y_sim, t(exp_i,(1:l)+time_index), 'makima' );

        end

        clear y2b dy2 HA2

        for kk = 1:l
             y2b(:,kk)      =  mean(y2(:,:,kk),2);
             HA2(:,:,kk)    =  y2(:,:,kk) - y2b(:,kk)*ones(1,q);
             HA2(:,:,kk)    =  HA2(:,:,kk)*TTinv;
             dy2(:,kk)      =  mean(yth(:,time_index+kk)) - y2b(:,kk);
        end

        E2(:,1:l)           =  squeeze(y2);
        aux2                =  zeros(q,l);

        for kk = 1:l
             [R(kk),~]      =  robustcov(yth(:,time_index+kk));
             Rinv           =  linsolve(R(kk),eye(size(R(kk))));
             aux1(:,:,kk)   =  (HA2(:,:,kk)'*Beta(kk)*Rinv*HA2(:,:,kk))/(q-1);
             aux2(:,kk)     =  aux2(:,kk) + (HA2(:,:,kk)'*Beta(kk)*Rinv*(dy2(:,kk)))/(q-1);

        end

        GGinv               =  eye(q) + sum(aux1,3);
        GG                  =  linsolve(GGinv,eye(size(GGinv)));
        b                   =  linsolve(GGinv,sum(aux2,2));
        dx1                 =  A10*b + A10*GG*pinv(A10'*A10)*A10'*(x10-x1);
        x1                  =  x1 + dx1;
        TT                  =  sqrtm(GG);

        x_est(:,jj)         =  x1;
        timestep_time(jj)   =  toc;
        E_est(:,:,jj)       =  E1;
        A1                  =  A10*TT;
        E1                  =  x1*ones(1,q) + A1;

        if min(E1(2*NT+NTM+8,:))<=0

            Idx = (E1(2*NT+NTM+8,:))<=0;
            E1(2*NT+NTM+8,Idx) = E1(2*NT+NTM+8,Idx)-min(E1(2*NT+NTM+8,:))+1e-6;

        end
            disp(['ended iteration ',num2str(jj),' with norm(dx1) = ', ...
                num2str(norm(dx1)), ' and norm(dy2) = ' num2str(norm(dy2)), ' at ',num2str(timestep_time(jj)),' seconds for the ' num2str(exp_i) 'th experiments '])

        jj = jj+1;
    end
    
    % Perform covaraince inflation:
    if j<j_max
        for j1 = 1:size(A1,1)
            if std(A1(j1,:)) ==0
                std_prior(j1) =0;
            else

                std_post(j1) = sqrt(robustcov(A1(j1,:)));
            end
        end

        RTPS              = 1+CI_theta*(std_prior-std_post)./std_post;
        RTPS(isnan(RTPS)) = 1;
        RTPS(isinf(RTPS)) = 1;

        for j1 = 1:size(A1,1)
            A1(j1,:) = mvnrnd(0,(RTPS(j1).*std_post(j1))^2,q);
        end

        Input_prior =  false;
        E1 = x1*ones(1,q) + A1;
    end

    x(:,:,j+1) = E1;


end


% Posterior ensembles:
G_post                 =   E1(2*NT+NTM+7,:)*P_inf;
mu_post                =   E1(2*NT+NTM+8,:)/Uc*(P_inf*R0);
alpha_post             =   (E1(2*NT+NTM+10,:));
[sigma_post,mu_post]   =   robustcov([G_post; mu_post; alpha_post]');
P_post.mu              =   mu_post;
P_post.sigma           =   sigma_post;


c = clock;
run_time = toc


end