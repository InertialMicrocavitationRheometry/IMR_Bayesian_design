function [mu_NMC] = IMR_EIG(Design, Model_all, N, sigma_all )


% addpath('/Users/tianyichu/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/GT/IMR/IMR_data_assimilation-main/IMR_DA')


N_model = size(Model_all,1);

Y = [];


for j = 1:N_model

    

    Model_j = Model_all{j};


    model = Model_j{1};

    P_prior = Model_j{3};

    model_prob = Model_j{4};

    disp([model ' model with ' num2str(round(100*model_prob),'%.2f') '% probability']);

    mu_theta=P_prior.mu; sigma_theta = P_prior.sigma;


    theta_params = mvnrnd(mu_theta,sigma_theta,round(N*model_prob));
    
    if model_prob>1e-2

    [t(j,:),y2,U2,dt] = IMR_simulations(theta_params,Design,model,60);


        % Y(:,:,j) =  [ y2 U2*dt];

    % Y(:,:,j) =  y2 ;

    Y(end+1:end+round(N*model_prob),:) = y2;

    end

end

t_exp = t(1,:);


% Y =  reshape(  permute(Y,[1 3 2]), size(Y,1)*size(Y,3),size(Y,2));


% data = {theta_params,Design,y2};  
% save('sim_data.mat','data','-v7.3')


%%


if ~isempty(sigma_all)

    
    Design_observed = zeros(size(sigma_all,2),2);
    sigma_observed  = zeros(size(sigma_all,2),81);

    for j1 = 1:size(sigma_all,2)

        Design_observed(j1,:) = sigma_all{j1}{1};
        sigma_observed(j1,:)  = sigma_all{j1}{2};
    end

    sigma_w = zeros(1,81);

    for k1 = 1:81
        % sigma_w(k1) = (interp2(Design_observed(:,1),Design_observed(:,2),sqrt(sigma_observed(:,k1)), Design(:,1),Design(:,2),  'makima' )).^2;

        F_interp = scatteredInterpolant(Design_observed(:,1),Design_observed(:,2),sqrt(sigma_observed(:,k1)),'natural','linear');

        sigma_w(k1) = (F_interp(Design(:,1), Design(:,2)));
    end
 
else

   sigma_w = linspace(1e-6,0.025,81);

end


 
% [mu_NMC_j,~,] = NMC_est(Y, [], diag([sigma_w.^2  sigma_w.^2 ]),'data');


 [mu_NMC_j,~,] = NMC_est(Y, [], diag(sigma_w(1:61).^2 ),'data');


mu_NMC = mean(mu_NMC_j);



end