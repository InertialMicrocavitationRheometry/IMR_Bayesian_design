function [x_opt, y_opt] = BayOpts_IMR(Model_all,obs,xrange,N,sigma_w)

%  IMR-based Bayesian optimal experimental design (BOED)

%  Input:   Model_all  --  All the constutitve models and their corresponding material properties under consideration
%           obs        --  Number of BO trials
%           xrange     --  searching range for the optimal design
%           N          --  sample size for approximating the EIG
%           sigma_w    --  variance of the synthetic noise
            f = @(x) IMR_EIG(x, Model_all, N, sigma_w );           % Target function: EIG

%  Output:  x_opt      --  the optimal design
%           y_opt      --  the optimal EIG

%% Bayesian Optimization 
clc;

[Xgrid,Ygrid]    =  meshgrid(linspace(xrange(1,1),xrange(1,2)),...
    linspace(xrange(2,1),xrange(2,2)));
XYgrid           =  [Xgrid(:), Ygrid(:)];
XEI_all          =  zeros(obs+1,size(Xgrid,1),size(Xgrid,2));
std_all          =  zeros(obs+1,size(Xgrid,1),size(Xgrid,2));
y_pred_all       =  zeros(obs+1,size(Xgrid,1),size(Xgrid,2));

%%  Load exsisting data or initialize with 10 random observations

save_name  = 'results_2models.mat';
load_data  = false;

if load_data
    load(save_name);                                             % load data from exisiting evaluations
    x      = data{1};
    y_true = data{2};
else
    Nstart = 10;                                                 % initial samples                         
    x      = [rand(Nstart,1)*diff(xrange(1,:))+xrange(1,1),...
              rand(Nstart,1)*diff(xrange(2,:))+xrange(2,1)];     % Random x positions
    y_true = zeros(size(x,1),1);
    for t = 1:size(x,1), y_true(t) = f(x(t,:)); end              % Evaluate the function for the initial samples
    data   = {x, y_true, y_pred_all,std_all,XEI_all};
    save(save_name,'data','-v7.3')
end

%%  Optimization process

for j = 0:obs

    disp(['Design parameters: Rmax = ' num2str(x(end,1)) ', Req = ' num2str(x(end, 2)) ', EIG = ' num2str(y_true(end))])

    % Gaussian regression,  ardmatern52 kernel was recommended in https://proceedings.neurips.cc/paper/2012/file/05311655a15b75fab86956663e1819cd-Paper.pdf

    if j == 0
        mdl = fitrgp(x,y_true,'Sigma',0.01,'ConstantSigma',true,...
            'KernelFunction','ardmatern52');
    else
        mdl = fitrgp(x,y_true,'KernelFunction','ardmatern52',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus','Very_true_max_idxse',0));
    end

    XYgrid              =  [Xgrid(:), Ygrid(:)];
    [y_pred,std]        =  predict(mdl,XYgrid);
    y_pred_all(j+1,:,:) =  reshape(y_pred,size(Xgrid));

    %% Expected Improvement

    xi                  =  0.01;  % Exploration-exploitation parameter (High xi = more exploration, Low xi = more exploitation)
    deviation           =  y_pred - max(y_true) - xi; 
    EI                  =  (std ~= 0).*(std.*normpdf(deviation./std)+deviation.*normcdf(deviation./std));
    [~,EI_idx]          =  max(EI); 
    x(end+1,:)          =  XYgrid(EI_idx,:);           
    [y_true(end+1)]     =  f(x(end,:));   
    XEI_all(j+1,:,:)    =  reshape(EI,size(Xgrid));
    std_all(j+1,:,:)    =  reshape(std,size(Xgrid));

end

[y_pred_max,y_pred_max_idx] = max(y_pred);
[y_true_max,y_true_max_idx] = max(y_true); 

fprintf('Bayesian Optimization\n');
fprintf('%s (estimated):\n\ty(%.6f,%.6f) = %.6f\n','Maximum',XYgrid(y_pred_max_idx,1),XYgrid(y_pred_max_idx,2),y_pred_max);
fprintf('%s (observed):\n\ty(%.6f,%.6f) = %.6f\n', 'Maximum',x(y_true_max_idx,1),x(y_true_max_idx,2),y_true_max);

x_opt   = x(y_true_max_idx,:);
y_opt   = y_true(y_true_max_idx);

end


