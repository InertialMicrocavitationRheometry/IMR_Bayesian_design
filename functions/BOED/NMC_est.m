function [mu_NMC_j,log_P_LKH] = NMC_est(F,theta,sigma_w,varargin)
 
% 
% nested Monte-Carlo estimator for the expected information gain


% Inputs: forward model, q = F(theta,d) + w
%         Prior distribution of theta, p(theta) or mu_theta and sigma_theta
%         design configuration, d
%         measurements y with dimension Nx * M;
%         distribution of the noise, p(w) or sigma_w
%         range for the model parameters theta_range = [theta_min thta_max];

%%  

%  Sample theta from the Prior distribution p(theta)


%  Calculate the results from the forward model



if nargin == 4 
    F_type  = varargin{1};
else
    F_type  = 'model';
end


if strcmpi(F_type,'model')
    
    q_model = F(theta);
    
elseif strcmpi(F_type,'data')
    
    q_model = F;
    
end



% size(theta)
% size(q_model)

[N,Nx] = size(q_model);



% 
% norm_fun     = @(mu, Sigma, x)  (1/(sqrt((2*pi)^(length(x))*det(Sigma))))+(-((x-mu)/(2*Sigma))*(x-mu)' );

 % norm_fun_log = @(mu, Sigma, x)  log(1/(sqrt((2*pi)^(Nx)*det(Sigma))))+(-((x-mu)/(2*Sigma))*(x-mu)' );

 norm_fun_log = @(mu, Sigma, x)  +(-((x-mu)/(2*Sigma))*(x-mu)' );


%%  draw samples using the likelihood function


y = zeros(N,Nx);

for n = 1:N

    y(n,:) = mvnrnd(q_model(n,:),sigma_w,1); 

end


%% Calculate the likelihood



% compute the all the likelihood functions

log_P_LKH   = zeros(N,N);


for j1 = 1:N


    for j2 = 1:N


    log_P_LKH(j1,j2) =   (norm_fun_log(q_model(j2,:), sigma_w,y(j1,:)))';

    end

     % log_P_prior(j1,1) =   (norm_fun_log(P_prior.mu, P_prior.sigma, theta(j1,:)))';

end




log_P_EVD = log(mean(exp(log_P_LKH),2));



%%
% mu_NMC_j = zeros(1,N);


 % log_P_post = log_P_prior + diag(log_P_LKH) -log_P_EVD;


mu_NMC_j   = diag(log_P_LKH) -log_P_EVD;



mu_NMC = mean(mu_NMC_j);




end
