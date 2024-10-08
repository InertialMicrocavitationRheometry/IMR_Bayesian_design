addpath('functions/BOED')
addpath('functions/DA')
addpath('functions/IMR simulation')

%% Set up parallel enviroment

poolobj = parpool('Processes', 24);
fprintf('Number of workers: %g\n', poolobj.NumWorkers);

%% Set up an underlying model

model_true = 'fung';   % Quadratic-law Kelvin--Voight
theta_true = [2770 0.186 0.48];
save_name  = 'results_design.mat';

%%  Initialize the prior distribution

%%%%%%%%%%%%%%%%%%%%%%---Model 1: NHKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_prior_1.mu      =  [15090 0.209 0];
P_prior_1.sigma   =  ([4350 0 0; 0 0.18 0; 0 0 0]).^2;
model_prob_1      =  0.5;
Model_1_prior     =  {'NeoHook', 'G, mu, alpha', P_prior_1, model_prob_1};

%%%%%%%%%%%%%%%%%%%%%%%---Model 2: qKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_prior_2.mu      =  [2770 0.286 0.28];
P_prior_2.sigma   =  ([300 0 0; 0 0.186 0; 0 0 0.48]).^2;
model_prob_2      =  0.5;
Model_2_prior     =  {'fung', 'G_inf, mu, alpha', P_prior_2, model_prob_2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Model_design      =  {Model_1_prior; Model_2_prior};
count             =  1;
model_prob_all(count,:) = [model_prob_1, model_prob_2];
Model_all{count}  = {Model_design; [];[]};

while count<30
    %%  Bayesian optimization

    N            =  1000;                     % EIG sample size
    obs          =  15;                       % BO evaluation numbers
    xrange       =  [100 1000; 0.14 0.3];     % Optimization range

    if count <4
        [Design_opt, EIG_opt] = BayOpts_IMR(Model_design,obs,xrange,N,[]);
    else
        [Design_opt, EIG_opt] = BayOpts_IMR(Model_design,obs,xrange,N,sigma_all);
    end

   disp(['Design #' num2str(count) ', Optimal design at We = ' num2str(round(Design_opt(1)),'%i') ', Req = ' num2str(Design_opt(2),'%.2f')])

    %% generate the data for the optimal design

    data_type     =  'sim'; % 'sim' or 'exp'
    data_set      =  'simulation_data';
    data_filepath =  'simulation_data/';
    data_filename =  ['sim_data', '_We',num2str(round(Design_opt(1)),'%i'),'_Req', ...
    num2str(Design_opt(2),'%.2f') '.mat']; % name of file containing R vs T data
    save_info     =  {data_type, data_set, data_filepath, data_filename};
    
    % default = 0.05, perform IMR simulations up to t = 4.

    [t,yth,Uth,dt] = IMR_simulations(theta_true,Design_opt,model_true,80,'auto',save_info);   
     
    % save error signatures
    sigma_w           = zeros(1,81);
    for kk = 1:80
        sigma_w(kk+1) = robustcov(yth(:,1+kk));
    end
    sigma_w(1)        = 1e-6;
    sigma_all{count}  = {Design_opt, sigma_w};

    %%  Data assimilation
    clear q_model

    N_model        =  size(Model_design,1);

    for j = 1:N_model
        Model_j    =  Model_design{j};
        if j == N_model
            % fixed range of quasistatic shear modulus for qKV
            Model_j{3}.mu(1)         = 2770;
            Model_j{3}.sigma(1,1)    = 300^2;
        end
        [x,E2,Model_post{j},params]  = IMR_DA( Model_j, save_info);
        q_model(j,:,:)               = squeeze(E2);
    end

    %% Model prbability

    Q_model                          =  {squeeze(q_model(1,:,:)); squeeze(q_model(2,:,:))};
    q_data                           =  yth(:,(1:size(q_model,3))+1);
    [model_prob_all(count+1,:)]      =  Model_prob_est(Q_model, q_data);

    %%  update the Prior

    % %%%%%%%%%%%%%%%%%%%%%%%---Model 1: NHKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Model_1_post                     =  Model_post{1};
    P_post_1                         =  Model_1_post{3};
    P_prior_1.mu(1:2)                =  ((P_prior_1.mu(1:2)/P_prior_1.sigma(1:2,1:2))+(P_post_1.mu(1:2)/P_post_1.sigma(1:2,1:2)))/((inv(P_prior_1.sigma(1:2,1:2))+inv(P_post_1.sigma(1:2,1:2))) ) ;
    P_prior_1.sigma(1:2,1:2)         =  inv(inv(P_prior_1.sigma(1:2,1:2))+inv(P_post_1.sigma(1:2,1:2)));
    model_prob                       =  mean(model_prob_all(1:count+1,1));
    Model_1_prior                    =  {'NeoHook', 'G, mu, alpha', P_prior_1, model_prob};

    disp(['Model: NeoHook, Posterior mean = ' num2str(P_prior_1.mu,'%.3f')])
    disp(['Model: NeoHook, Posterior std = ' num2str(sqrt(diag(P_prior_1.sigma))','%.3f')])

    %%%%%%%%%%%%%%%%%%%%%%%%%---Model 2: qKV----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Model_2_post                     =  Model_post{2};
    P_post_2                         =  Model_2_post{3};
    P_prior_2.mu                     =  ((P_prior_2.mu/P_prior_2.sigma)+(P_post_2.mu/P_post_2.sigma))/((inv(P_prior_2.sigma)+inv(P_post_2.sigma)) ) ;
    P_prior_2.sigma                  =  inv(inv(P_prior_2.sigma)+inv(P_post_2.sigma));
    model_prob                       =  mean(model_prob_all(1:count+1,2));
    Model_2_prior                    =  {'fung', 'G_inf, mu, alpha', P_prior_2, model_prob};

    disp(['Model: qKV, Posterior mean = ' num2str(P_prior_2.mu,'%.3f')])
    disp(['Model: qKV, Posterior std = ' num2str( sqrt(diag(P_prior_2.sigma))','%.3f')])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Model_design                     =  {Model_1_prior; Model_2_prior};
    Model_all{count+1}               =  {Model_design; Design_opt; EIG_opt};

    save(save_name,'Model_all','-v7.3')

    count = count +1;

end



