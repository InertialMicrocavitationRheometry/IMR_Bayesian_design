% This script determines initial X0 and x0
%{
global tspan R0 NT NTM Pext_type Pext_Amp_Freq Tgrad Cgrad model G G1 ...
    mu t0 neoHook nhzen sls linkv k chi fom foh We Br A_star B_star ...
    Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star De deltaY ...
    yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm x0_true N
%}
%***************************************
% Extract viscoelastic parameters from stuct
G = visco_params.G;
G1 = visco_params.G1;
mu = visco_params.mu;
alpha = visco_params.alpha;
lambda_nu = visco_params.lambda_nu;

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
tspan_star = tspan/t0;
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

tau_del = [];

x0_true = [X0,Br,foh,G,mu,De,alpha,lambda_nu,est_params];
N = length(x0_true);