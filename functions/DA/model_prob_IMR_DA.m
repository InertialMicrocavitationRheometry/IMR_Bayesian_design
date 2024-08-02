%% Model probability test for IMR_DA





%%   load models

addpath ./results


for j1 = 1:11

load('DA_En4D_NeoHook_exp.mat');

q_model_1 = squeeze(E2(:,:,[20]));

q_model_1 = reshape( permute(q_model_1,[1,3,2]),[],size(q_model_1,2));

load('DA_En4D_fung_exp.mat');


q_model_2 = squeeze(E2(:,:,[20]));

q_model_2 = reshape( permute(q_model_2,[1,3,2]),[],size(q_model_2,2));


Q_model = {q_model_1(:,:); q_model_2(:,:)};





%%   load exp data

% load('sim_data_noisy.mat');
% 
% 
% % q_data = mvnrnd(y, (0.01*eye(Nx)).^2+diag( [ linspace(1e-4,3e-2,Nx/2) linspace(1e-4,3e-2,Nx/2)  ]).^2 ,N) ;
% 
% 
% q_data = yth(:,(1:size(q_model_2,2))+1);


data_type = 'exp'; % 'sim' or 'exp'
data_set = 'SoftPA_nobeads';
data_filepath = (['data/PA_05%_0.03%/']);
data_filename = 'PA_05%_0.03%.mat'; % name of file containing R vs T data


N_exp = 20;
import_data_exp_PA;

q_data = yth([1:6, 7:20],(1:size(q_model_2,2))+1);


 [model_prob_all(j1, :)] = Model_prob_est(Q_model, q_data)

end


figure;

subplot(1,2,1)

 plot((q_model_1(:,:))','r:'); hold on;
plot((q_data(:, (1:size(q_model_2,2)) ))','kx'); axis square; hold on

title('Model 1: NeoHook')
xlabel('t^*')
ylabel('R^*')


subplot(1,2,2)


plot((q_model_2(:,:))','b:'); hold on;
plot((q_data(:,(1:size(q_model_2,2))))','kx'); axis square; hold on

title('Model 2: Gen qKV')
xlabel('t^*')
ylabel('R^*')





%%


figure; 

subplot(2,1,1)

plot(model_prob_all(1,:),'r'); hold on; plot(model_prob_all(2,:),'b');


legend('Model 1: Neohook','Model 2: qKV')

ylabel('log_{10}(Marginal Likelihood)')
xlabel('Number of measurements')


subplot(2,1,2)

plot(model_prob_all(1,:)-model_prob_all(2,:),'m');

ylabel('\Delta')
xlabel('Number of measurements')



%%  


figure;

plot(model_prob_all(1:11,1)-model_prob_all(1:11,2),'m');

ylabel('log(Bayes factor): NeoHook/Gen qKV')
xlabel('Number of measurements')