function [t_exp,y2,U2,dt] = IMR_simulations(theta_params,Design,model,time_window,varargin)

%%  Running IMR for a set of prior paramters

%  Input:   theta_params  --  Distributions of material properties
%           Design        --  Design parameters (e.g. Rmax, Req)
%           model         --  Constitutive models:  'neoHook','nhzen','sls','linkv','fung','fung2','fungexp','fungnlvis'
%           time_window   --  size of simulation window
%           varagin       --  Options for saving data and adding noises


tic


%% Initialize and import data


q         =  size(theta_params,1);
G1        =  1e9;                
lambda_nu =  .1;


%% Modeling parameters



NT = 240; % Amount of nodes inside the bubble
NTM = 240; % Amount of nodes outside the bubble
Pext_type = 'IC'; % Type of external forcing
ST = 0.056; % (N/m) Liquid Surface Tension

Tgrad = 1; % Thermal effects inside bubble
Tmgrad = 1; % Thermal effects outside bubble
Cgrad = 1; % Vapor diffusion effects
comp = 1; % Activates the effect of compressibility (0=Rayleigh-Plesset, 1=Keller-Miksis)
disp_timesteps = 1; % displays timesteps in En4D run (for debugging)

% following should not be changed (untested):
disptime = 0; % 1 = display simulation time
Dim = 0; % 1 = output variables in dimensional form


est_params = [];


%%


S          =  0.056; %0.056; % JY!!! % 0.056; % (N/m) Liquid Surface Tension
P_inf      =  101325; % (Pa) Atmospheric Pressure
rho        =  1060; % (Kg/m^3) Liquid Density
R0         =  Design(1)/P_inf*(2*S);
Req        =  Design(2);


opts.POSDEF=  true;
opts.SYM   =  true;
l          =  time_window; % size of data assimilation window
time_index =  1;


% Pre-allocating memory:

n          =  80;
tspan_star =  4;
t          =  linspace(0,tspan_star,n+1);
Uc         =  sqrt(P_inf/rho);
t0         =  R0/Uc;
t1         =  t(time_index);
t2         =  t(time_index+l);
t_exp      =  t(time_index+(0:l));
dt         =  (t(2)-t(1));
vars       =  cell(q,1);
tau_del    =  cell(q,1);



parfor memb = 1:q

    G     = theta_params(memb,1); 
    mu    = theta_params(memb,2);  
    alpha = theta_params(memb,3);

    if mu<=0
        mu = 1e-6;
    end

    %%


    addpath ./IMR-vanilla/functions
    addpath ./IMR-vanilla/src


    %% Guess for parameters

    %NT = 500; %400 % Ammount of nodes inside the bubble (>=500 is a good to start)
    %NTM = 500; %20 % Ammount of nodes outside the bubble (>=10 is a good to start)
    %Pext_type = 'IC'; % Type of external forcing. Refer to RP_Cav

    % Find Req and calculate initial partial pressure
    % Needed parameters
    %ST = 0.056; % (N/m) Liquid Surface Tension
    Pmt_temp = IMRcall_parameters(R0,G,G1,mu); % Calls parameters script
    Ca = Pmt_temp(5); Re = Pmt_temp(6);
    P_inf = Pmt_temp(19); T_inf = Pmt_temp(20);

    % Req = mean(yth(end-5:end)); % take mean of end of sim to be Req
    P_guess = (P_inf + (2*ST)/(Req*R0) - Pvsat(T_inf))*(Req^3);

    Pext_Amp_Freq =[P_guess 0]; % [ Pressure ; Freq ], Freq = 0 for IC
    %Pext_Amp_Freq = [100 0];

    % Simulation parameters
    %disptime = 0; % 1 = display simulation time
    %Tgrad = 1; % Thermal effects inside bubble

    %if exist('Tmgrad')
    %    % don't modify if already set
    %else
    %    Tmgrad = 1;% Thermal effects outside bubble
    %end

    %Cgrad = 1;  % Vapor diffusion effects
    %Dim = 0; % Output variables in dimensional form
    %comp = 1; % Activates the effect of compressibility

    %
    %G1 = inf;
    %G = G_guess;
    %mu = mu_guess;

    %% Determine initial state vector based on parameters

    % Load Parameters :
    Pmt = IMRcall_parameters(R0,G,G1,mu); % Calls parameters script
    k = Pmt(1); chi = Pmt(2); fom = Pmt(3); foh = Pmt(4); Ca = Pmt(5);
    Re = Pmt(6); We = Pmt(7); Br = Pmt(8); A_star = Pmt(9); B_star = Pmt(10);
    Rv_star = Pmt(11); Ra_star = Pmt(12); P0_star = Pmt(13); t0 = Pmt(14);
    C0 = Pmt(15); L = Pmt(16); L_heat_star = Pmt(17); Km_star = Pmt(18);
    P_inf = Pmt(19); T_inf = Pmt(20); C_star = Pmt(21); De = Pmt(22); rho = Pmt(23);

    %****************************************

    % Material Choice
    neoHook = 0;
    nhzen = 0;
    sls = 0;
    linkv = 0;
    fung = 0; fung2 = 0; fungexp = 0; fungnlvis = 0;
    if strcmp(model,'neoHook') == 1
        neoHook = 1;
    elseif strcmp(model,'fung') == 1
        fung = 1;
    elseif strcmp(model,'fung2') == 1
        fung2 = 1;
    elseif strcmp(model,'fungexp') == 1
        fungexp = 1;
    elseif strcmp(model,'fungnlvis') == 1
        fungnlvis = 1;
    elseif strcmp(model,'nhzen') == 1
        nhzen = 1;
    elseif strcmp(model,'sls') == 1
        sls = 1;
    elseif strcmp(model,'linkv') == 1
        linkv = 1;
    else
        nhzen = 1;
    end


    % Needed to account for fast diffusion
    P0_star = P0_star - (1-Cgrad)*Pvsat(1*T_inf)/P_inf;

    % When we assume water vapor undergoes infinitely fast mass diffusion
    % the vapor pressure is constant and P is the pressure of
    % non-condesible gas

    %******************************************
    % Creates finite difference matrices
    D_Matrix_T_C = Finite_diff_mat(NT,1,0);
    DD_Matrix_T_C = Finite_diff_mat(NT,2,0);
    D_Matrix_Tm = Finite_diff_mat(NTM,1,1);
    DD_Matrix_Tm = Finite_diff_mat(NTM,2,1);
    %******************************************

    %******************************************
    % Create spatial nodes

    % Inside the bubble
    N = NT-1;
    deltaY = 1/N;
    i = 1:1:N+1;
    yk = ((i-1)*deltaY)';

    % Outside the bubble
    Nm = NTM-1;
    deltaYm = -2/Nm;
    j = 1:1:Nm+1;
    xk = (1+(j-1)*deltaYm)';
    yk2 = ((2./(xk+1)-1)*L+1);

    %******************************************
    % Initial Conditions
    R0_star = 1;
    U0_star = 0;  % Change as needed
    %Z10 = 0;
    S0 = 0;
    Tau0 = zeros(1,NT);
    C0 = C0*ones(1,NT);
    Tm0 = ones(1,NTM);
    if strcmp(Pext_type,'ga')
        dt_star = Pext_Amp_Freq(2)/t0;
        w_star = Pext_Amp_Freq(3)*t0;
    end

    % Need to modify intial conditions for the Out-of-Equilibrium Rayleigh
    % Collpase:
    if strcmp(Pext_type,'IC')
        Pv = Pvsat(1*T_inf)/P_inf;
        P0_star = Pext_Amp_Freq(1)/P_inf + Cgrad*Pvsat(1*T_inf)/P_inf;
        % Need to recalculate intital concentration
        theta = Rv_star/Ra_star*(P0_star-Pv)/Pv; % mass air / mass vapor
        C0 = 1/(1+theta);

        % Calculate the equilibrium radii ratio for initial stress state:
        if Req == 0
            [REq,~,~] = IMRCalc_Req(R0, Tgrad, Cgrad, Pext_Amp_Freq(1), G, G1, mu);
        end
        REq = Req;
        %REq = 1; %removed 6/15/16 by Jon
        C0 = C0*ones(1,NT);
        %U0_star = -1*(1-P0_star)/(C_star); %Intitial velocity due to shockwave
        U0_star = 0;

        if sls == 1 || linkv == 1
            S0 = -4/(3*Ca)*(1-REq^3);
        elseif nhzen == 1 || neoHook == 1
            S0 = -1/(2*Ca)*(5-REq^4-4*REq);
        elseif fung == 1 || fungnlvis == 1
            S0 = -(1-3*alpha)*(5 - 4*REq - REq^4)/(2*Ca) - ...
                2*alpha*(-27/40 - 1/8*REq^8 - 1/5*REq^5 -1*REq^2 + 2/REq)/(Ca);
        elseif fung2 == 1
            S0 = 2*(1-3*alpha+4.5*alpha^2)*(-5/4 + 1*REq + 1/4*REq^4)/(Ca) + ...
                2*(27/40*alpha-221/90*alpha^2 + alpha^2/24*REq^12 + alpha^2/18*REq^9 + (alpha-3*alpha^2)/8*REq^8 + ... )
                2*alpha^2/6*REq^6 + (alpha-3*alpha^2)/5*REq^5 + 2*alpha^2/3*REq^3 + (2*alpha-6*alpha^2)/2*REq^2 - ...
                2*alpha^2*log(1/REq) - (2*alpha-6*alpha^2)/REq - 2/3*alpha^2/REq^3)/(Ca);
        elseif fungexp == 1
            tempbeta = linspace(1/REq,1,1e5);
            tempS = 2*(tempbeta.^-5 + tempbeta.^-2) .* exp(alpha*(tempbeta.^-4+2*tempbeta.^2-3))/(Ca);
            % figure, plot(tempbeta,tempS);
            S0 = (trapz(tempbeta,tempS));
        end

    end

    X0 = [R0_star U0_star P0_star S0 Tau0 C0 Tm0];

    % tau_del = [];

    x0_true = [X0,Br,foh,G,mu,De,alpha,lambda_nu,est_params];

    %%
    % xi(:,memb) = x0_true;

    vars{memb} = {NT Pext_type Pext_Amp_Freq disptime Tgrad Tmgrad ...
        Cgrad comp t0 neoHook nhzen sls linkv k chi fom foh We Br A_star ...
        B_star Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star ...
        De deltaY yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
        D_Matrix_Tm DD_Matrix_Tm tspan_star NTM rho R0 fung fung2 fungexp fungnlvis};


    [t_memb{memb}, EE{memb},~] = f_IMR(t1,t2,x0_true,vars{memb},tau_del{memb});


    t_sim = t_memb{memb};
    y_sim = EE{memb}(:,1);
    U_sim = EE{memb}(:,2);

    y2(memb,:) = interp1(t_sim,y_sim, t_exp, 'makima' );
    U2(memb,:) = interp1(t_sim,U_sim, t_exp, 'makima' );


end



c = clock;
run_time = toc


%%  save data (add syhthetic noise if necessary)

if nargin >4 && ~isempty(varargin{1})

    add_noise = true;
    sigma_w   = varargin{1};
    if strcmp(varargin{1},'auto')
        sigma_w = 0.02*abs(mean(y2,1)-1) + linspace(0,0.025,81);
    else
        sigma_w = interp1(linspace(0,4,81),sqrt(sigma_w), t_exp, 'makima' );
    end

    q = 100;
    y2 =  mvnrnd(y2, diag(sigma_w.^2) ,q) ;
    U2 =  mvnrnd(y2, diag(sigma_w.^2) ,q) ;
else
    add_noise = false;
    sigma_w   = 0;
end



if nargin >5
    save_info  = varargin{2};
    save_data = true;
    data_type = save_info{1}; % 'sim' or 'exp'
    data_set  = save_info{2};
    data_filepath = save_info{3};
    data_filename = save_info{4}; % name of file containing R vs T data

else
    save_data = false;
end




if save_data

    for j = 1:q

        sim_data(j).R0     =  R0;
        sim_data(j).t0     =  t0;
        sim_data(j).Roft   =  y2(j,:)*R0;
        sim_data(j).R_norm =  y2(j,:);
        sim_data(j).U_norm =  U2(j,:);
        sim_data(j).Req    =  Req;
        sim_data(j).model  =  model;
        sim_data(j).Design =  Design;

        if add_noise
            sim_data(j).t      =  t_exp*t0;
            sim_data(j).t_norm =  t_exp;
            sim_data(j).theta  =  theta_params(1,:);
        else
            sim_data(j).t      =  t_memb{j}*t0;
            sim_data(j).t_norm =  t_memb{j};
            sim_data(j).theta  =  theta_params(j,:);
        end
    end


save([ data_filepath,data_filename],'sim_data','-v7.3')


end


end





