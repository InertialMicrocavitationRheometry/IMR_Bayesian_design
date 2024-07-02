function [x_opt, y_opt] = BayOpts_2D_IMR(Model_all,obs,xrange,N,sigma_w)

%% Bayesian Optimization for 2D opti problem

% close all; clc;

clc;

angle = [45 45];
% Nstart = 1;    % Initial no. of observations


f = @(x) IMR_EIG(x, Model_all, N, sigma_w );


tobemax = true;

[Xfine,Yfine] = meshgrid(linspace(xrange(1,1),xrange(1,2)),...
    linspace(xrange(2,1),xrange(2,2)));


xyfine = [Xfine(:), Yfine(:)];


zmax = 10;


XEI_all = zeros(obs+1,size(Xfine,1), size(Xfine,2));
sd_all = zeros(obs+1,size(Xfine,1),size(Xfine,2));
y_pred_all = zeros(obs+1,size(Xfine,1),size(Xfine,2));

save_name = 'results_2models.mat';


load_data = false;

if load_data

    load(save_name)

    x = data{1};
    y_true = data{2};

else

    Nstart = 10;
    x = [rand(Nstart,1)*diff(xrange(1,:))+xrange(1,1),...
        rand(Nstart,1)*diff(xrange(2,:))+xrange(2,1)]; % Random x positions
    y_true = zeros(size(x,1),1);



    for t = 1:size(x,1), y_true(t) = f(x(t,:)); end     % Sample the func at x

    data = {x, y_true, y_pred_all,sd_all,XEI_all};

    save(save_name,'data','-v7.3')

end



for j = 0:obs


    disp(['Design parameters: Rmax = ' num2str(x(end,1)) ', Req = ' num2str(x(end, 2)) ', EIG = ' num2str(y_true(end))])


if j == 0
        mdl = fitrgp(x,y_true,'Sigma',0.01,'ConstantSigma',true,...
        'KernelFunction','ardmatern52');
else
    
  mdl = fitrgp(x,y_true,'KernelFunction','ardmatern52',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus','Verbose',0));
end

close all;

    % ardmatern52 kernel was recommended in https://arxiv.org/pdf/1206.2944.pdf

    xyfine = [Xfine(:), Yfine(:)];
    [y_pred,sd] = predict(mdl,xyfine);

    y_pred_all(j+1,:,:) = reshape(y_pred,size(Xfine));



     figure;   

    subplot(222);
    surf(Xfine,Yfine,reshape(y_pred,size(Xfine))); 	% GPR prediction
    shading interp; hold on; view(angle);
    scatter3(x(:,1),x(:,2),y_true,10,'g','filled',...
        'MarkerEdgeColor','k');                     % Plot seen data
    title(sprintf('No. of observed points: %d',length(x)));
    axis([xrange([1 3 2 4]), zlim]); box on; hold off;
    drawnow;

    subplot(223);
    surf(Xfine,Yfine,reshape(sd,size(Xfine))); 	% Uncertainty (std. dev.)
    shading interp; hold on; view(angle);
    axis([xrange([1 3 2 4]), 0, max(2,max(sd))]);
    title('Uncertainty'); box on; hold off;
        xlabel('Rmax')
    ylabel('Req')


    %% Expected Improvement
    % This EI is from http://krasserm.github.io/2018/03/21/bayesian-optimization/

    xi = 0.01;  % Exploration-exploitation parameter (greek letter, xi)
    % High xi = more exploration
    % Low xi = more exploitation (can be < 0)

    if tobemax, d = y_pred - max(y_true) - xi; % (y - f*) if maximization
    else,       d = min(y_true) - y_pred - xi; % (f* - y) if minimiziation
    end

    EI = (sd ~= 0).*(d.*normcdf(d./sd) + sd.*normpdf(d./sd));

    [eimax,posEI] = max(EI); xEI = xyfine(posEI,:);
    x(end+1,:) = xEI;               %#ok<SAGROW> Save xEI as next


    [y_true(end+1)] = f(x(end,:));    %#ok<SAGROW> Sample the obj. at xEI
    

    XEI_all(j+1,:,:) = reshape(EI,size(Xfine));
    sd_all(j+1,:,:) = reshape(sd,size(Xfine));

data = {x, y_true, y_pred_all,sd_all,XEI_all};

save(save_name,'data','-v7.3')





    subplot(224);
    surf(Xfine,Yfine,reshape(EI,size(Xfine)));
    shading interp; hold on;
    plot3(xEI([1 1]),xrange(2,:),eimax*[1 1],'--k','LineWidth',1.5);
    plot3(xrange(1,:),xEI([2 2]),eimax*[1 1],'--k','LineWidth',1.5);
    title('Expected Improvement'); grid on; hold off;
    view(angle); axis(xrange([1 3 2 4])); box on;
    xlabel('Rmax')
    ylabel('Req')



    subplot(222); hold on; xlabel('x1'); ylabel('x2');
    scatter3(x(end,1),x(end,2),zmax,'m','filled','MarkerEdgeColor','k');
    plot3(xEI([1 1]),xrange(2,:),zmax*[1 1],'--k','LineWidth',1.5);
    plot3(xrange(1,:),xEI([2 2]),zmax*[1 1],'--k','LineWidth',1.5);
    hold off; axis(xrange([1 3 2 4])); view(angle); pause(0.5);
    xlabel('Rmax')
    ylabel('Req')

end

 [ae,be] = max(y_pred);
 [ao,bo] = max(y_true); str = 'Maximum';


fprintf('Bayesian Optimization\n');
fprintf('  %s (estimated):\n\ty(%.6f,%.6f) = %.6f\n',...
    str,xyfine(be,1),xyfine(be,2),ae);
fprintf('  %s (observed):\n\ty(%.6f,%.6f) = %.6f\n',...
    str,x(bo,1),x(bo,2),ao);

x_opt = x(bo,:);

y_opt = y_true(bo);

% P_post_exp = P_post_exp{bo};

% P_post_exp = [];

 end

% data = {x, y_true, y_pred_all,sd_all,XEI_all};
% 
% save('results.mat','data','-v7.3')

