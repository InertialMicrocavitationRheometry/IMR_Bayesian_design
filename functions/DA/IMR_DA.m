function [x,E2,Model_params_post,params]= IMR_DA( Model_params, load_info)

%%  model =  'neoHook','nhzen','sls','linkv','fung','fung2','fungexp','fungnlvis'
%   q = 48; % Number of ensemble members


% clear all

%% Data import
% see import_data_exp.m file for details on data loading. The following is
% meant for data in the same format as RofT data from Selda's experiments.
% If the format is different, the import_data_exp.m file will need to be
% modified accordingly

        data_type = load_info{1}; % 'sim' or 'exp'
        data_set  = load_info{2};
        data_filepath = load_info{3};
        data_filename = load_info{4}; % name of file containing R vs T data


%dataset = 2; % data set number (if needed - old data format)

num_peaks = 2; % number of 'peaks' to assimilate in radius data
               % (peaks = collapse points as in Estrada paper)

import_data;



%% Run main for corresponding method:
   



  [x,E2,P_post,params]  = main_En4D_peaks_par(Model_params,t,yth,R0_all,Req_all,tspan_all,peak_time_idx);



  % E1 = squeeze(x);



%%


 Model_params_post{1} =  Model_params{1};
 Model_params_post{2} =  Model_params{2};



 Model_params_post{3} = P_post; 



end
%%