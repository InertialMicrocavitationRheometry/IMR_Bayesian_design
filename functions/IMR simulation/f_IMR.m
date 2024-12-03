function [t,X,tau_del] = f_IMR(ti_star,tf_star,xi,vars,tau_del)

global NT Pext_type Pext_Amp_Freq disptime Tgrad Tmgrad ...
    comp t0 neoHook nhzen sls linkv k chi fom foh We Br A_star ...
    B_star Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star ...
    De deltaY yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm Ca Re

NT= vars{1};Pext_type=vars{2};Pext_Amp_Freq=vars{3};disptime=vars{4};
Tgrad=vars{5};Tmgrad=vars{6};Cgrad=vars{7};comp=vars{8};t0=vars{9};
neoHook=vars{10};nhzen=vars{11};sls=vars{12};linkv=vars{13};k=vars{14};
chi=vars{15};fom=vars{16};foh=vars{17};We=vars{18};Br=vars{19};
A_star=vars{20};B_star=vars{21};Rv_star=vars{22};Ra_star=vars{23};
L=vars{24};L_heat_star=vars{25};Km_star=vars{26};P_inf=vars{27};
T_inf=vars{28};C_star=vars{29};De=vars{30};deltaY=vars{31};yk=vars{32};
deltaYm=vars{33};xk=vars{34};yk2=vars{35};Pv=vars{36};REq=vars{37};
D_Matrix_T_C=vars{38};DD_Matrix_T_C=vars{39};D_Matrix_Tm=vars{40};
DD_Matrix_Tm=vars{41}; tspan_star = vars{42}; NTM = vars{43}; rho = vars{44};
R0 = vars{45}; fung = vars{46}; fung2 = vars{47}; fungexp = vars{48};
fungnlvis = vars{49};

%************************************************
% March equations in time

Uc = sqrt(P_inf/rho);

Br = xi(2*NT+NTM+5);
foh = xi(2*NT+NTM+6);

Ca = P_inf./(xi(2*NT+NTM+7));
Re = (P_inf*R0)./((xi(2*NT+NTM+8)).*Uc);

% Ca = 1./(xi(2*NT+NTM+7));
% Re = 1./((xi(2*NT+NTM+8)));


De = xi(2*NT+NTM+9);
alpha = (xi(2*NT+NTM+10));
lambda_nu = xi(2*NT+NTM+11);


% restrictions on dynamics to keep quantities physical:

xi(3) = (xi(3));
xi(5+NT:4+(2*NT)) = max(xi(5+NT:4+(2*NT)),0);

options = odeset('RelTol',1e-5,'AbsTol',1e-8);
%[~ ,X] = ode23tb(@bubble, [ti_star tf_star], xi(1:end-2)',options);

 [t ,X] = ode23tb(@bubble, [ti_star tf_star], xi(1:2*NT+NTM+4)',options);

% [t ,X] = ode23tb(@bubble, [ti_star tf_star], xi(1:2*NT+NTM+4)');

% xf = [X(end,:)';Br;foh;xi(2*NT+NTM+7:end)];
% 
% xf(3) = (xf(3));


%%
%*************************************************************************
% Nested function; ODE Solver calls to march governing equations in time
% This function has acess to all parameters above

    function dxdt = bubble(t,x)
        % Break x vector into indv. values
        R = x(1); % Bubble wall Radius
        U = x(2); % Bubble wall velocity
        P = x(3); % Internal pressure
        %P = max(0.1,x(3)); % artificial fix to avoid crash
        %Z1 = x(4); % Stress integral
        S = x(4);
        Tau = x(5:(NT+4));
        %Tau = max(-0.7,x(5:(NT+4))); % artificial fix to avoid crash
        C = x((NT+5):(2*NT+4));
        Tm = x((2*NT+5):end);
        
        if (disptime == 1)
            disp(t/tspan_star);
        end
        
        %*********Solves for boundary condition at the wall**************
        if (Tmgrad == 1)
            if t/tspan_star> 0.001
                %Might need to tune 0.001 for convergence:


                guess= real(-.001 + tau_del(end)); %+tau_del(end);

                if isnan(guess) || isinf(guess)
                    guess
                end

                prelim  = fzero(@(x) Boundary(x,Tm,Tau,C,P),guess);
            else
                guess = -.0001;
                prelim  = fzero(@(x) Boundary(x,Tm,Tau,C,P),guess);
            end
        else
            prelim = 0;
        end
        
        %****************************************************************
        % Sets value at boundary conditions
        tau_del = [tau_del prelim]; % Commented out here for simplification
        
        %tau_del_end = prelim;
        
        Tau(end) = prelim;
        T = TW(Tau);
        Tm(1) = T(end);
        % TW(prelim)
        
        % Calculated variables
        K_star = A_star*T+B_star;
        C(end) =  CW(T(end),P);
        
        Rmix = C*Rv_star + (1-C)*Ra_star;
        
        % Gets variables that are not directly calculated as outputs
        
        % JY!!! ************ Update values because of changes of Cp *******
        % Also update Cp-> Br, Dm -> Foh
        % Already called: rho = 1060; % (Kg/m^3) Liquid Density
        Rnondim = P_inf/(rho*T_inf);
        Cp = Rmix(end)*k/(k-1);
        Uc = sqrt(P_inf/rho);
        Br = Uc^2./(Cp*T_inf);
        A = 5.3e-5; % (W/m-K^2)Thermal Conductivity coeff
        B = 1.17e-2; %(W/m-K)Thermal Conductivity coeff
        K_infy = A*T_inf+B;
        Dm = Km_star*(K_infy) ./ (rho*Cp);
        foh = Dm/(Uc*R0);
        
        % This is commented out in this version for simplification
        
        %Tdel = [Tdel T(end)];
        %tdel = [tdel t];
        %Cdel = [Cdel C(end)];
        
        
        % The following was commented out at first:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Set external pressure
        %if (Pext_type == 'sn')
        %    Pext =  -Pext_Amp_Freq(1)/P_inf*sin( Pext_Amp_Freq(2)*t*t0) ;
        %    P_ext_prime = -Pext_Amp_Freq(2)*t0*Pext_Amp_Freq(1)/P_inf...
        %        *cos( Pext_Amp_Freq(2)*t*t0) ;
        %elseif (Pext_type == 'RC')
        %    Pext = Pext_Amp_Freq(1)/P_inf ;
        %    P_ext_prime = 0;
        %elseif (Pext_type == 'RG')
        %    Pext = -Pext_Amp_Freq(1)/P_inf ;
        %    P_ext_prime = 0;
        %elseif (Pext_type == 'ip')
        %    Pext = -Pext_Amp_Freq(1)/P_inf*...
        %        (1-heaviside(t-Pext_Amp_Freq(2)/t0)) ;
        %    P_ext_prime = 0;
        %end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if strcmp(Pext_type,'IC')
            Pext = 0;
            P_ext_prime = 0;
        elseif strcmp(Pext_type , 'ga')
            Pext = pf(t)/P_inf;
            P_ext_prime = pfdot(t)/P_inf;
        end
        
        % *****************************************
        % Create derivative terms
        
        % Temp. field of the gas inside the bubble
        DTau  = D_Matrix_T_C*Tau;
        DDTau = DD_Matrix_T_C*Tau;
        
        % Concentration of vapor inside the bubble
        DC  = D_Matrix_T_C*C;
        DDC = DD_Matrix_T_C*C;
        
        % Temp. field in the material
        DTm = D_Matrix_Tm*Tm;
        DDTm = DD_Matrix_Tm*Tm;
        
        %***************************************
        % Internal pressure equation
        %pdot = 3/R*(Tgrad*chi*(k-1)*DTau(end)/R - k*P*U +...
        %    + Cgrad*k*P*fom*Rv_star*DC(end)/( T(end)*R* Rmix(end)* (1-C(end)) ) );
        pdot = 3/R*(Tgrad*chi*(k-1)*DTau(end)/R - k*P*U +...
            + Cgrad*k*P*fom*Rv_star*DC(end)/( R* Rmix(end)* (1-C(end)) ))  ;
        % *****************************************
        
        %***************************************
        % Temperature of the gas inside the bubble
        %U_vel = (chi/R*(k-1).*DTau-yk*R*pdot/3)/(k*P);
        %first_term = (DDTau.*chi./R^2+pdot).*(K_star.*T/P*(k-1)/k);
        %second_term = -DTau.*(U_vel-yk*U)./R;
        
        %Tau_prime = first_term+second_term;
        %Tau_prime(end) = 0;
        %Tau_prime = Tau_prime*Tgrad;
        
        % Temperature of the gas inside the bubble
        % U_vel = (chi/R*(k-1).*DTau-yk*R*pdot/3)/(k*P);
        U_vel = (chi/R*(k-1).*DTau-yk*R*pdot/3)/(k*P) + fom/R*(Rv_star-Ra_star)./Rmix.*DC;
        first_term = ((DDTau ).*chi./R^2+pdot).*(K_star.*T/P*(k-1)/k);
        second_term = -DTau.*(U_vel-yk*U)./R;
        third_term = fom./(R.^2) *(Rv_star-Ra_star)./Rmix .* DC .*DTau;
        
        Tau_prime = first_term + second_term + third_term;
        if Tmgrad == 0
            Tau_prime(end) = 0;
        else
            %    Tau_prime(end) = K; %JY???
            Tau_prime(end) = 0; % What is K?
        end
        Tau_prime = Tau_prime*Tgrad;
        
        % *****************************************
        
        %***************************************
        % Vapor concentration equation
        %U_mix = U_vel + fom/R*((Rv_star - Ra_star)./Rmix).*DC;
        %one = DDC;
        %two = DC.*(DTau./(K_star.*T)+((Rv_star - Ra_star)./Rmix).*DC );
        %three =  (U_mix-U.*yk)/R.*DC;
        %
        % C_prime = fom/R^2*(one - two) - three;
        % C_prime(end) = 0;
        % C_prime = C_prime*Cgrad;
        
        % Vapor concentration equation
        U_mix = U_vel; % + fom/R*((Rv_star - Ra_star)./Rmix).*DC;
        one = DDC;
        % % JY!!!  %
        %two = DC.*( -((Rv_star - Ra_star)./Rmix).*DC - DTau./(K_star.*T) );
        two = DC.*( -((Rv_star - Ra_star)./Rmix).*DC - DTau./sqrt(1+2*A_star.*Tau)./T );
        %tempdata = sqrt(1+2*A_star.*Tau);
        %disp('======')
        %[Tau(1:40),tempdata(1:40),K_star(1:40)]
        %disp('======')
        
        three = (U_mix-U.*yk)/R.*DC;
        
        % % JY!!! % C_prime = fom/R^2*(one - two) - three;
        C_prime = fom/R^2*(one+two) - three;
        C_prime(end) = 0;
        C_prime = C_prime*Cgrad;
        %*****************************************
        
        %***************************************
        % Material temperature equations
        %first_term = (1+xk).^2./(L*R).*(U./yk2.^2.*(1-yk2.^3)/2+foh/R.*((xk+1)/(2*L)-1./yk2)).* DTm;
        %second_term = foh/R^2.*(xk+1).^4/L^2.*DDTm/4;
        %third_term =  3*Br./yk2.^6.*(4/(3*Ca).*(1-1/R^3)+4.*U/(Re.*R)).*U./R;
        %Tm_prime = first_term+second_term+third_term;
        %Tm_prime(end) = 0; % Sets boundary condition on temp
        %Tm_prime(1) = 0; % Previously calculated;
        %Tm_prime = Tm_prime*Tmgrad; %Tmgrad makes this quantity zero
        
        % Material temperature equations
        first_term = (1+xk).^2./(L*R).*(U./yk2.^2.*(1-yk2.^3)/2+foh/R.*((xk+1)/(2*L)-1./yk2)).* DTm;
        second_term = foh/(R^2).*(xk+1).^4/L^2.*(DDTm)/4; %JY???
        % % JY!!! third_term =  3*Br./yk2.^6.*(4/(3*Ca).*(1-1/R^3)+4.*U/(Re.*R)).*U./R;
        
        third_term = 3*Br./yk2.^6.*(4.*U/(Re.*R)).*U./R;
        
        Tm_prime = first_term+second_term+third_term;
        Tm_prime(end) = 0; % Sets boundary condition on temp
        Tm_prime(1) = 0; % Previously calculated;
        Tm_prime = Tm_prime*Tmgrad; %Tmgrad makes this quantity zero
        %*****************************************
        
        %***************************************
        %{
        % Elastic stress in the material
        if linkv == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                S = -4/(3*Ca)*(1 - 1/Rst^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R/Rst^3 + 4/Re*U^2/R^2;
            else
                S = -4/(3*Ca)*(1 - 1/R^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R^4 + 4/Re*U^2/R^2;
            end
        elseif neoHook == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                S = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) - 4/Re*U/R ;
                Sdot =  -2*U/R*(1/Rst + 1/Rst^4)/Ca + 4/Re*U^2/R^2;
            else
                S = -(5 -4/R - 1/R^4)/(2*Ca) - 4/Re*U/R;
                Sdot =  -2*U*(1/R^2 + 1/R^5)/Ca + 4/Re*U^2/R^2;
            end
        elseif sls == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                Sdot = -S/De - 4*(1-1/Rst^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            else
                Sdot = -S/De - 4*(1-1/R^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            end
        elseif nhzen == 1
            if  strcmp(Pext_type, 'IC')
                Rst = R/REq;
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                if isinf(Sdot)
                    Rst=Rst+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                end
            else
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                if isinf(Sdot)||isnan(Sdot)
                    R = R+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                end
            end
        end
        %}
        
        % Elastic stress in the material
        if linkv == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                S = -4/(3*Ca)*(1 - 1/Rst^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R/Rst^3 + 4/Re*U^2/R^2;
            else
                S = -4/(3*Ca)*(1 - 1/R^3) - 4/Re*U/R;
                Sdot =  -4/Ca*U/R^4 + 4/Re*U^2/R^2;
            end
        elseif neoHook == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                % ====== Old for neoHook ======
                S = -(5 - 4/Rst - 1/Rst^4)/(2*Ca) - 4/Re*U/R ;
                Sdot =  -2*U/R*(1/Rst + 1/Rst^4)/Ca + 4/Re*U^2/R^2;
                % JY!!! "- 4/Re*udot/R" is added finally: Sdot = Sdot - SdotA*udot/R; % JY!!! Pay attention to here!
            else
                S = -(5 -4/R - 1/R^4)/(2*Ca) - 4/Re*U/R;
                Sdot =  -2*U*(1/R^2 + 1/R^5)/Ca + 4/Re*U^2/R^2;
            end
            
        elseif fung == 1 || fungnlvis == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                % ====== JY!!! for Fung model ======
                % alpha = 0.1;
                % ******** JY!!! First order Fung model approx ======
                S = -(1-3*alpha)*(5 - 4/Rst - 1/Rst^4)/(2*Ca) - ...
                    2*alpha*(-27/40 - 1/8/Rst^8 - 1/5/Rst^5 -1/Rst^2 + 2*Rst)/(Ca) - 4/Re*U/R;
                Sdot = -2*U/R*(1-3*alpha)*(1/Rst + 1/Rst^4)/Ca - ...
                    2*alpha*U/R*(1/Rst^8 + 1/Rst^5 + 2/Rst^2 + 2*Rst)/(Ca) + 4/Re*U^2/R^2;
                %JY!!! "- 4/Re*udot/R" is added later: Sdot = Sdot - SdotA*udot/R; % JY!!! Pay attention to here!
            else
                disp('Not finished in non-IC case for Fung model!')
            end
            % ====== JY!!! First order Fung G + first order mu model approx ======
            % % alpha = 0; lambda_nu = 0.001;
            % Lv = 1;
            %
            % zeta = linspace(-1,0.99,200);
            %
            % tempr = R*( (2./(1-zeta)-1)*Lv + 1 );
            % tempr0 = (tempr.^3+R0^3-R^3).^(1/3);
            % gammadot = -0.5*( 2*(tempr0.^2)./(tempr.^3) + 1./tempr0 ) *R^2./(tempr.^2) * U;
            % % figure; plot(zeta,gammadot); pause;
            % % tempmu = 1/Re .* heaviside(-abs(gammadot)+1/lambda_nu) .* (1-lambda_nu^2*(gammadot.^2));
            % tempmu = 1/Re .* exp(-lambda_nu^2*(gammadot.^2));
            % tempS = -12*2*tempmu*U/R .* (1-zeta).^2 ./ (2+(1-zeta)*(1/Lv-1)).^4 /(Lv^3);
            %
            % S = -(1-3*alpha)*(5 - 4/Rst - 1/Rst^4)/(2*Ca) - ...
            %     2*alpha*(-27/40 - 1/8/Rst^8 - 1/5/Rst^5 -1/Rst^2 + 2*Rst)/(Ca) + ...
            %     trapz(zeta,tempS);
            %
            % Sdot = -2*U/R*(1-3*alpha)*(1/Rst + 1/Rst^4)/Ca - ...
            %      2*alpha*U/R*(1/Rst^8 + 1/Rst^5 + 2/Rst^2 + 2*Rst)/(Ca) + 4/Re*U^2/R^2;
        elseif fung2 == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                % ******** JY!!! Second order Fung model approx ******
                S = 2*(1-3*alpha+4.5*alpha^2)*(-5/4 + 1/Rst + 1/4/Rst^4)/(Ca) + ...
                    2*(27/40*alpha-221/90*alpha^2 + alpha^2/24/Rst^12 + alpha^2/18/Rst^9 + (alpha-3*alpha^2)/8/Rst^8 + ...
                    2*alpha^2/6/Rst^6 + (alpha-3*alpha^2)/5/Rst^5 + 2*alpha^2/3/Rst^3 + (2*alpha-6*alpha^2)/2/Rst^2 - ...
                    2*alpha^2*log(Rst) - (2*alpha-6*alpha^2)*Rst - 2/3*alpha^2*Rst^3)/(Ca) - ...
                    4/Re*U/R;
                Sdot = 2*U/R*(1-3*alpha+4.5*alpha^2)*(-1/Rst-1/Rst^4)/(Ca) + ...
                    2*U/R*(-alpha^2/2/Rst^12 - alpha^2/2/Rst^9 - (alpha-3*alpha^2)/Rst^8 ...
                    -2*alpha^2/Rst^6 - (alpha-3*alpha^2)/Rst^5 - 2*alpha^2/Rst^3 - (2*alpha-6*alpha^2)/Rst^2 ...
                    -2*alpha^2 - (2*alpha-6*alpha^2)*Rst - 2*alpha^2*Rst^3)/(Ca) + 4/Re*U^2/R^2;
                %JY!!! "- 4/Re*udot/R" is added later: Sdot = Sdot - SdotA*udot/R;
            else
                disp('Not finished in non-IC case for Fung2 model!')
            end
            
        elseif fungexp == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                % ******** JY!!! The following original Fung model codes don't work. ********
                tempbeta = linspace(Rst,1,1e5);
                tempS = 2*(tempbeta.^-5 + tempbeta.^-2) .* exp(alpha*(tempbeta.^-4+2*tempbeta.^2-3))/(Ca);
                S = (trapz(tempbeta,tempS)) - 4/Re*U/R;
                Sdot = -2*U/REq*(Rst^-5 + Rst^-2)*exp(alpha*(Rst^(-4)+2*Rst^2-3))/(Ca) + ...
                    4/Re*U^2/R^2;
            else
                disp('Not finished in non-IC case for Fung2 model!')
            end
            
            
        elseif sls == 1
            if  strcmp(Pext_type,'IC')
                Rst = R/REq;
                Sdot = -S/De - 4*(1-1/Rst^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            else
                Sdot = -S/De - 4*(1-1/R^3)/(3*Ca*De) - 4/(Re*De)*U/R - 4*U/(Ca*R);
            end
            
        elseif nhzen == 1
            if  strcmp(Pext_type, 'IC')
                Rst = R/REq;
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                if isinf(Sdot)
                    Rst=Rst+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/Rst^4-4/Rst)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(Rst^3-1)*(3/14*Rst^3+Rst^2-3/(2*Rst)+2/(7*Rst^4));
                end
            else
                Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                    -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                if isinf(Sdot)||isnan(Sdot)
                    R = R+eps;
                    Sdot = -S/De - 1/(2*Ca*De)*(5-1/R^4-4/R)-4*U/(R*Re*De)...
                        -4*U/(R*Ca)/(R^3-1)*(3/14*R^3+R^2-3/(2*R)+2/(7*R^4));
                end
            end
        end
        
        %****************************************************
        
        % Equations of motion
          rdot = U;
        
        %         if (Tgrad == 0)
        %             P = P0_star*(1/R)^(3*k);
        %             pdot = -3*k*U/R*P;
        %         end
        
        Pv = (Pvsat(T(end)*T_inf)/P_inf);
        if comp == 0
            %Rayleigh-Plesset equation
            udot = (P + abs(1-Cgrad)*Pv  - 1 - Pext + S -1/(We*R) -1.5*U^2)/R;
            SdotA = 0;
        else
            % Keller-Miksis equation
            if linkv==1 || neoHook==1 || fung==1 || fung2==1 || fungexp==1 || fungnlvis==1
                SdotA = 4/Re;
            elseif sls==1 || nhzen==1
                SdotA = 0;
            end
            %if fungexp==0
            udot = ((1+U/C_star)...
                *(P  + abs(1-Cgrad)*Pv -1/(We*R) + S - 1 - Pext)  ...
                + R/C_star*(pdot+ U/(We*R^2) + Sdot -P_ext_prime ) ...
                - 1.5*(1-U/(3*C_star))*U^2)/((1-U/C_star)*R);%+JdotA/(C_star));
            
            %elseif fungexp==1
            %    part1 = Sdot;
            %    part2 = ((1+U/C_star)...
            %    *(P  + abs(1-Cgrad)*Pv -1/(We*R) + S - 1 - Pext)  ...
            %    + R/C_star*(pdot+ U/(We*R^2) + 0 -P_ext_prime ) ...
            %    - 1.5*(1-U/(3*C_star))*U^2)/((1-U/C_star)*R);
            %    Sdot = (1+4/Re/C_star)^(-1) * (part1-4/Re/R*part2);
            %    udot = (1+4/Re/C_star)^(-1) * (part2+R/C_star*part1);
            %end
        end
        % ****************************************
        if  fungnlvis==1
            Sdot = Sdot - SdotA/lambda_nu * ((3.4641016)^(lambda_nu-1)) * ((abs(udot)/R)^(lambda_nu));
        else
            Sdot = Sdot - SdotA*udot/R; % JY!!! Pay attention to here!
        end
        
        dxdt = [rdot; udot; pdot; Sdot; Tau_prime; C_prime; Tm_prime];
        if isreal(rdot)==0
            pause;
        end
    end
%*************************************************************************

% Other nested functions used to carry out repetetive calculations
% throughout the code
    function Tw = TW(Tauw)
        %calculates the temperature at the bubble wall as a fuction of \tau
        Tw = (A_star -1 + sqrt(1+2*Tauw*A_star)) / A_star;
    end

    function Cw = CW(Tw,P)
        % Calculates the concentration at the bubble wall
        
        %Function of P and temp at the wall
        theta = Rv_star/Ra_star*(P./(Pvsat(Tw*T_inf)/P_inf) -1);
        Cw = 1./(1+theta);
    end

    function Tauw = Boundary(prelim,Tm,Tau,C,P)
        % Solves temperature boundary conditions at the bubble wall
        % Create finite diff. coeffs.
        % Coefficients in terms of forward difference
        
        %    %Second order
        coeff = [-3/2 , 2 ,-1/2 ];
        Tm_trans = Tm(2:3);
        T_trans = flipud(Tau(end-2:end-1));
        C_trans = flipud(C(end-2:end-1));
        
        %   Could implement any order... sixth order example is shown below
        
        %     Sixth order
        %     coeff= [-49/20 ,6	,-15/2	,20/3	,-15/4	,6/5	,-1/6]; %Sixth order coeff
        %     Tm_trans = Tm(2:7);
        %     T_trans = flipud(Tau(end-6:end-1));
        %     C_trans = flipud(C(end-6:end-1));
        
        Tauw =chi*(2*Km_star/L*(coeff*[TW(prelim); Tm_trans] )/deltaYm) +...
            chi*(-coeff*[prelim ;T_trans] )/deltaY + Cgrad*...
            fom*L_heat_star*P*( (CW(TW(prelim),P)*(Rv_star-Ra_star)+Ra_star))^-1 *...
            (TW(prelim) * (1-CW(TW(prelim),P))  ).^(-1).*...
            (-coeff*[CW(TW(prelim),P); C_trans] )/deltaY;
        %************************************************************************
    end

% Gaussian pressure functions
% acoustic pressure
    function p = pf(t)
        if t<(dt_star-5*tw_star) || t>(dt_star+5*tw_star)
            p=0;
        else
            p = -Pext_Amp_Freq(1)*exp(-(t-dt_star).^2/tw_star^2);
        end
    end

% time derivative of acoustic pressure
    function pdot = pfdot(t)
        if t<(dt_star-5*tw_star) || t>(dt_star+5*tw_star)
            pdot=0;
        else
            pdot = 2*(t-dt_star)/tw_star^2*Pext_Amp_Freq(1).*exp(-(t-dt_star).^2/tw_star^2);
        end
        
    end

end
