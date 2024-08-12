% This file sets up an initial ensemble for experimental data setups

% State vector used: [R,U,P,S,Tau,C,Tm, 1/Ca, 1/Re]

addpath ./IMR-vanilla/functions
addpath ./IMR-vanilla/src

% Shuffle rng to ensure randomness of results
rng('shuffle')


%% Guess for parameters

%NT = 500; %400 % Ammount of nodes inside the bubble (>=500 is a good to start)
%NTM = 500; %20 % Ammount of nodes outside the bubble (>=10 is a good to start)
%Pext_type = 'IC'; % Type of external forcing. Refer to RP_Cav

% Find Req and calculate initial partial pressure
% Needed parameters
%ST = 0.056; % (N/m) Liquid Surface Tension
Pmt_temp = IMRcall_parameters(R0,G_guess,G1_guess,mu_guess); % Calls parameters script
Ca = Pmt_temp(5); Re = Pmt_temp(6);
P_inf = Pmt_temp(19); T_inf = Pmt_temp(20);

% Req = mean(yth(end-5:end)); % take mean of end of sim to be Req
P_guess = (P_inf + (2*ST)/(Req*R0) - Pvsat(T_inf))*(Req^3);

Pext_Amp_Freq =[P_guess 0]; % [ Pressure ; Freq ], Freq = 0 for IC
%Pext_Amp_Freq = [100 0];



%% Determine initial state vector based on parameters
initialize

%% Create initial ensemble


% Custom spread:
x_init = x0_true';

if Input_prior

   Caspread    = 0;
   Respread    = 0;
   alphaspread = 0;

end

spread = [Rspread; Uspread; Pspread; Sspread; ones(NT,1)*tauspread; ...
    ones(NT,1)*Cspread; ones(NTM,1)*Tmspread; Brspread; fohspread; ...
    Caspread; Respread; Despread; alphaspread; lambda_nuspread];

xi = (1 + spread .* randn(N,q)) .* repmat(x_init,1,q) + ...
    repmat([0;0;0;0;zeros(2*NT+NTM,1);0;0;0;0;0;0;0],1,q) .* randn(N,q);



%%  using truncated distribution 

if Input_prior

    xi(2*NT+NTM+7,:)      =   G_prior';
    xi(2*NT+NTM+8,:)      =   mu_prior';
    xi(2*NT+NTM+10,:)     =   alpha_prior';

else

    pd_G                  =   makedist('Normal','mu',G_guess,'sigma',G_guess*Caspread);
    t_G                   =   truncate(pd_G,1e-10,inf);                                  % truncated Gaussian distribution

    if Caspread ==0
        xi(2*NT+NTM+7,:)  =   G_guess;
    else
        xi(2*NT+NTM+7,:)  =   random(t_G,[1, q]);
    end

    pd_mu                 =   makedist('Normal','mu',mu_guess,'sigma',mu_guess*Respread);
    t_mu                  =   truncate(pd_mu,1e-10,inf);
    xi(2*NT+NTM+8,:)      =   random(t_mu,[1, q]);
    pd_alpha              =   makedist('Normal','mu',alpha_guess,'sigma',alpha_guess*alphaspread);

    if alphaspread ==0  
        xi(2*NT+NTM+10,:) = alpha_guess;
    else
        xi(2*NT+NTM+10,:) = random(t_alpha,[1, q]);
    end

end

%%

xi(3,:) = log(xi(3,:));
xi(3,:) = log(x0_true(3));     % constrain the pressure to be positive
Uc      = sqrt(P_inf/rho);
xi      = [xi(1:2*NT+NTM+6,:);...
           (xi(2*NT+NTM+7,:))./P_inf; ...
           ((xi(2*NT+NTM+8,:)).*Uc)./(P_inf*R0);...
           xi(end-2,:);...
           (xi(end-1,:));...
           xi(end,:)]; % for now

