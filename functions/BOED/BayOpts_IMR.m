function [x_opt, y_opt] = BayOpts_IMR(Model_all,obs,xrange,N,sigma_w)

%  IMR-based Bayesian optimal experimental design (BOED)

%  Input:   Model_all  --  All the constutitve models and their corresponding material properties under consideration
%           obs        --  Number of BO trials
%           xrange     --  searching range for the optimal design
%           N          --  sample size for approximating the EIG
%           sigma_w    --  variance of the synthetic noise
            f = @(x) IMR_EIG(x, Model_all, N, sigma_w );           % Target function: EIG

%  Output:  x_opt      --  the optimal design
%           y_opt      --  associated optimal EIG

%% Bayesian Optimization 
clc;

[Xfine,Yfine]    =  meshgrid(linspace(xrange(1,1),xrange(1,2)),...
    linspace(xrange(2,1),xrange(2,2)));
xyfine           =  [Xfine(:), Yfine(:)];
XEI_all          =  zeros(obs+1,size(Xfine,1),size(Xfine,2));
sd_all           =  zeros(obs+1,size(Xfine,1),size(Xfine,2));
y_pred_all       =  zeros(obs+1,size(Xfine,1),size(Xfine,2));

%%  Load exsisting data or initialize with 10 random observations

save_name = 'results_2models.mat';

load_data = false;

if load_data
    load(save_name);                                     % load data from exisiting evaluations
    x      = data{1};
    y_true = data{2};
else

    Nstart = 10;                                         % initial samples                         
    x = [rand(Nstart,1)*diff(xrange(1,:))+xrange(1,1),...
        rand(Nstart,1)*diff(xrange(2,:))+xrange(2,1)];   % Random x positions
    y_true = zeros(size(x,1),1);

    for t = 1:size(x,1), y_true(t) = f(x(t,:)); end      % Sample the func at x
    data = {x, y_true, y_pred_all,sd_all,XEI_all};
    save(save_name,'data','-v7.3')

end

%%  Optimization process

for j = 0:obs

    disp(['Design parameters: Rmax = ' num2str(x(end,1)) ', Req = ' num2str(x(end, 2)) ', EIG = ' num2str(y_true(end))])

    % Gaussian regression

    if j == 0
        mdl = fitrgp(x,y_true,'Sigma',0.01,'ConstantSigma',true,...
            'KernelFunction','ardmatern52');
    else
        mdl = fitrgp(x,y_true,'KernelFunction','ardmatern52',...
            'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus','Verbose',0));
    end

    % ardmatern52 kernel was recommended in https://arxiv.org/pdf/1206.2944.pdf

    xyfine              =  [Xfine(:), Yfine(:)];
    [y_pred,sd]         =  predict(mdl,xyfine);
    y_pred_all(j+1,:,:) =  reshape(y_pred,size(Xfine));

    %% Expected Improvement

    xi               =  0.01;  % Exploration-exploitation parameter (High xi = more exploration, Low xi = more exploitation)
    d                =  y_pred - max(y_true) - xi; 
    EI               =  (sd ~= 0).*(d.*normcdf(d./sd) + sd.*normpdf(d./sd));
    [~,posEI]        =  max(EI); 
    xEI              =  xyfine(posEI,:);
    x(end+1,:)       =  xEI;               
    [y_true(end+1)]  =  f(x(end,:));   
    XEI_all(j+1,:,:) =  reshape(EI,size(Xfine));
    sd_all(j+1,:,:)  =  reshape(sd,size(Xfine));

end

[ae,be] = max(y_pred);
[ao,bo] = max(y_true); str = 'Maximum';

fprintf('Bayesian Optimization\n');
fprintf('  %s (estimated):\n\ty(%.6f,%.6f) = %.6f\n',...
    str,xyfine(be,1),xyfine(be,2),ae);
fprintf('  %s (observed):\n\ty(%.6f,%.6f) = %.6f\n',...
    str,x(bo,1),x(bo,2),ao);

x_opt   = x(bo,:);
y_opt   = y_true(bo);

end


