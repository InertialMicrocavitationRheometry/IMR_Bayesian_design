function [mu_NMC] = IMR_EIG(Design, Model_all, N, sigma_all )

%  IMR-based expected information gain (EIG) estimator

%  Input:   Design     --  A given design point
%           Model_all  --  All the constutitve models and their corresponding material properties under consideration
%           N          --  sample size for approximating the EIG
%           sigma_all  --  Avavilable variances of the synthetic noise at different designs

%  Output:  mu_NMC     --  Nested Monte-Carlo EIG estimator

%%  Collecting samples from all the avaiable models

N_model = size(Model_all,1);

Y = [];

for j = 1:N_model
    
    % Input model setups

    Model_j      =  Model_all{j};
    model        =  Model_j{1};
    P_prior      =  Model_j{3};
    model_prob   =  Model_j{4};
    mu_theta     =  P_prior.mu; 
    sigma_theta  =  P_prior.sigma;
    theta_params =  mvnrnd(mu_theta,sigma_theta,round(N*model_prob));

    disp([model ' model with ' num2str(round(100*model_prob),'%.2f') '% probability']);
    
    if model_prob>1e-2

        [t(j,:),y2,U2,dt] = IMR_simulations(theta_params,Design,model,60);
        Y(end+1:end+round(N*model_prob),:) = y2;

    end

end

t_exp = t(1,:);

%%  Obtain the variances from the known error signatures


if ~isempty(sigma_all)

    Design_observed = zeros(size(sigma_all,2),2);
    sigma_observed  = zeros(size(sigma_all,2),81);

    for j1 = 1:size(sigma_all,2)
        Design_observed(j1,:) = sigma_all{j1}{1};
        sigma_observed(j1,:)  = sigma_all{j1}{2};
    end

    sigma_w = zeros(1,81);
    for k1 = 1:81
        F_interp = scatteredInterpolant(Design_observed(:,1),Design_observed(:,2),sqrt(sigma_observed(:,k1)),'natural','linear');
        sigma_w(k1) = (F_interp(Design(:,1), Design(:,2)));
    end
 
else

   sigma_w = linspace(1e-6,0.025,81);

end

%%  Evalute the EIG

[mu_NMC_j,~,] =  NMC_est(Y, [], diag(sigma_w(1:61).^2 ),'data');
mu_NMC        =  mean(mu_NMC_j);

end

