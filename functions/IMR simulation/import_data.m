% This file exports data in vector yth, it must be re-written for each data
% set to ensure the data is formatted correctly

%dataset = 8; %between 1-20 for water case
%addpath ./IMR-vanilla/functions
%addpath ./IMR-vanilla/src
% load('/scratch/jspratt/EnKS/exp_data/11kPa_PA/RofTdata.mat')


%% new data format
exp_data0 = load([data_filepath,data_filename]);

exp_data = exp_data0.sim_data;

clear yth t R_exp t_exp noise 

N_exp = size(exp_data,2);

 for exp_i = 1:N_exp

    R_exp = exp_data(exp_i).Roft(1:1:end);
    t_exp(exp_i,:) = exp_data(exp_i).t_norm(1:1:end);

    R0_all(exp_i) = exp_data(exp_i).R0;

    % fixing shape if needed
    if length(R_exp(1,:)) == 1
        R_exp = R_exp';
    end
    if length(t_exp(1,:)) == 1
        t_exp = t_exp';
    end
    
   
    t_exp(exp_i,:) = t_exp(exp_i,:)- (t_exp(exp_i,1)*2-t_exp(exp_i,2));

    % [~,max_index] = max(R_exp);

    max_index = 1;

    yth(exp_i,:) = R_exp(max_index:max_index+80)./R0_all(exp_i);
    t(exp_i,:)   = t_exp(exp_i,max_index:max_index+80);


    % Req_all(exp_i) = mean(yth(end-5:end)); % take mean of end of sim to be Req


     Req_all(exp_i) = exp_data(exp_i).Req;

     t(exp_i,:) = t(exp_i,:)- (t(exp_i,1));

     % R0_all(j) = exp_data(j).R0;


    % delet NaNs
    kk = 1;
    for jj = 1:size(yth,2)
        if isnan(yth(exp_i,jj))
        else
            yth(exp_i,kk) = yth(exp_i,jj);
            t(exp_i,kk) = t(exp_i,jj);
            kk = kk + 1;
        end
    end
    yth(exp_i,:) = yth(exp_i,1:kk-1);
    t(exp_i,:)   = t(exp_i,1:kk-1);
    %

    tspan_all(exp_i) = t(exp_i,end)-t(exp_i,1);
    n = size(t,2)-1;

    if exist('l') == 0
        l = n;
    end
    timesteps = n+1;

    % Find peak_time
    peak_indices = find(islocalmin(mean(yth(:,:),1)));
    peak_time_idx(exp_i) = peak_indices(num_peaks);
    

 end

