function [model_prob_all] = Model_prob_est(Q_model, q_data,  varargin)

%
% Model probability


% Inputs: forward model, q = F(theta,d) + w
%         Prior distribution of theta, p(theta) or mu_theta and sigma_theta
%         design configuration, d
%         measurements y with dimension Nx * M;
%         distribution of the noise, p(w) or sigma_w
%         range for the model parameters theta_range = [theta_min thta_max];

%%



[N,Nx] = size(q_data);


N_model = size(Q_model,1);



norm_fun     = @(mu, Sigma, x)  (1/(sqrt((2*pi)^(Nx)*det(Sigma))))*exp(-((x-mu)/(2*Sigma))*(x-mu)' );

 % norm_fun_log = @(mu, Sigma, x)  log(1/(sqrt((2*pi)^(Nx)*det(Sigma))))+(-((x-mu)/(2*Sigma))*(x-mu)' );


 % norm_fun_log = @(mu, Sigma, x)   -1/2*(log(det(Sigma))+Nx*log(2*pi)+(((x-mu)/(Sigma))*(x-mu)' ) );

norm_fun_log = @(mu, Sigma, x)   -1/2*(+(((x-mu)/(Sigma))*(x-mu)' ) );


model_prob_all = zeros(N_model,1);




%%  Distribution of Data

sigma_data = zeros(Nx);
mu_data    = zeros(1,Nx);

for j = 1:Nx

[sigma_data(j,j),mu_data(j)] = robustcov(q_data(:,j));

end


    %    mu_data = mean(q_data,1);
    % 
    % 
    % sigma_data = (0.015*eye(Nx)).^2+diag((linspace(1e-4,1.77e-2,Nx)).^2);



% sigma_data = (0.01*eye(Nx)).^2+diag( [ linspace(1e-4,3e-2,Nx/2) linspace(1e-4,3e-2,Nx/2)  ]).^2;


     % 

     P_EVD = zeros(N,1);

    % P_EVD = 0;


    for j1 = 1:N


        P_EVD(j1) =  (norm_fun_log(q_data(j1,:), sigma_data, mu_data))';

                % P_EVD(j1) =  (norm_fun_log(y(1,:), sigma_data, q_data(j1,:)))';

    end

    % P_EVD = exp(P_EVD);

    % P_EVD =  (norm_fun(mu_data, sigma_data, mu_data))';

%% Calculate the likelihood


    % compute the all the likelihood functions

     % P_LKH_log_max = zeros(:,N_model);


for k = 1:N_model


       q_model = Q_model{k};

       [M,~] = size(q_model);


    P_LKH   = zeros(1,M);


   for j2 = 1:M

      for j1 = 1:N


         % P_LKH(j1,j2) =   (norm_fun(q_model(j2,:), sigma_data, q_data(j1,:)))';

          P_LKH(1,j2) =  P_LKH(1,j2) + 1/N*(norm_fun_log(q_data(j1,:) , sigma_data, q_model(j2,:) ))';

      end

   end

   % P_LKH = exp(P_LKH);


    % P_EVD = mean(P_LKH,2);


       % model_prob_all(k) = mean(mean((P_LKH)./P_EVD));


 model_prob_all(k) = log10(mean(exp(P_LKH)));


% P_LKH_log(k,:) = P_LKH/N;


end



  model_prob_all = (10.^model_prob_all)./sum(10.^model_prob_all);

end