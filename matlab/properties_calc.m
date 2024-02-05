clear
clc 
clf

%% Parameters for reading file
out = [];
len = 1001; % thermo size 
n_variables = 18; % number of properties at thermo 
part = zeros(n_variables, len);
ipart = 0;

% Constants
%kb = 1.987204259e-3; %boltzmann constant kcal/mol/k
kb = 1.38064852e-23;            % [J/K]
n_a = 6.02214076e23; %avogrado number
kcalmol2J = 4184/n_a;
A32m3 = 1e-30;

% Simulation parameters - change it accordingly necessity
n_atoms = 400;          %% Molecules number in simulation
%Ti = 200 : 10 :390;     %% Temperature Input
mole_mass = 78.11;       %% g/mol

%run numbers x thermo number plot
T = zeros(20,1001);
V = zeros(20,1001);
rho = zeros(20,1001);
U = zeros(20,1001);
Etotal = zeros(20,1001);
H = zeros(20,1001);
H2 = zeros(20,1001);
V2 = zeros(20,1001);
VH = zeros(20,1001);
UH = zeros(20,1001);

%% Read LAMMPS file and get variables

for i = 1:20
   filename = [num2str(i),'.log'];
    fid = fopen(filename,'r');
    if fid < 0, error('Cannot open file'); end

    while 1  % Infinite loop
      s = fgets(fid);
      if ischar(s)
        data = sscanf(s, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g', n_variables);
        if length(data) == n_variables
          ipart = ipart + 1;
          part(:, ipart) = data;
          if ipart == len
            out = cat(2, out, part);
            ipart = 0;
          end
        end
      else  % End of file:
        break;
      end
    end

    out = cat(2, out, part(:, 1:ipart));
    %Temperature
    T(i,:) = part(2,:);

    %rho
    rho(i,:) = part(4,:);
    
    %Volume
    V(i,:) = part(5,:)*A32m3;%/(n_atoms/n_a); %m3/mol
    
    %Total Energy
    %E(i,:) = part(7,:); %kcal/mol
    Etotal(i,:) = part(8,:); %kcal/mol energia total

    %Internal Energy
    U(i,:) = part(6,:); %kcal/mol energia conf

    %Enthalpy
    H(i,:) = (Etotal(i,:)*4184+101325*V(i,:))/4184*kcalmol2J; %Enthalpy: J/mol
    Hconf(i,:) = (U(i,:)*4184+101325*V(i,:))/4184*kcalmol2J; %Enthalpy: J/mol - pro caso OPLS usar Emol e nao PotEng

    %Others
    %U =  U*kcalmol2J;
    VH(i,:) = V(i,:).*H(i,:);
    VHconf(i,:) = V(i,:).*Hconf(i,:);
    V2(i,:) = V(i,:).^2;
    H2(i,:) = H(i,:).^2;
    Hconf2(i,:) = Hconf(i,:).^2;
    UH(i,:) = U(i,:).*H(i,:);
    UHconf(i,:) = U(i,:).*Hconf(i,:);
    
end

%% Mean and Std Calculate
mean_T = transpose(mean(transpose(T)));                            %K
mean_rho = transpose(mean(transpose(rho)));                        %angstrom^3
mean_V = transpose(mean(transpose(V)));                            %m3
mean_H = transpose(mean(transpose(H)));                            %Enthalpy: J/mol
mean_Hconf = transpose(mean(transpose(Hconf)));                            %Enthalpy: J/mol
mean_U = transpose(mean(transpose(U)))*kcalmol2J;                  %Internal Energy U: J/mol
mean_Etotal = transpose(mean(transpose(Etotal)))*kcalmol2J;        %Total Energy: J/mol
mean_VH = transpose(mean(transpose(VH)));
mean_VHconf = transpose(mean(transpose(VHconf)));
mean_V2 = transpose(mean(transpose(V2)));
mean_H2 = transpose(mean(transpose(H2)));
mean_Hconf2 = transpose(mean(transpose(Hconf2)));
mean_UH = transpose(mean(transpose(UH)))*kcalmol2J;
mean_UHconf = transpose(mean(transpose(UHconf)))*kcalmol2J;

std_V = transpose(std(transpose(V))); 
std_V2 = transpose(std(transpose(V2)));
std_H = transpose(std(transpose(H))); 
std_H2 = transpose(std(transpose(H2))); 
std_U= transpose(std(transpose(U))); 
std_UH= transpose(std(transpose(UH))); 
std_VH = transpose(std(transpose(VH))); 
% std_V_mbar = transpose(std(transpose(V))); 
% std_V2_mbar = transpose(std(transpose(V2)));
% std_H_mbar = transpose(std(transpose(H))); 
% std_H2_mbar = transpose(std(transpose(H2))); 
% std_U_mbar= transpose(std(transpose(U)))*kcalmol2J; 
% std_UH_mbar= transpose(std(transpose(UH)))*kcalmol2J; 
% std_VH_mbar = transpose(std(transpose(VH))); 

%% Read MBAR files and NIST file                                              

V  = load('resultV.txt')*A32m3;
U  = load('resultU.txt')*kcalmol2J;
H  = load('resultH.txt')*kcalmol2J;
VH = load('resultVH.txt')*A32m3*kcalmol2J;
UH = load('resultUH.txt')*kcalmol2J*kcalmol2J;
V2 = load('resultV2.txt')*A32m3*A32m3;
H2 = load('resultH2.txt')*kcalmol2J*kcalmol2J;
T_mbar  = load('T.txt');  
Ti = load ('temp.inp');
Ti = transpose(Ti);

% Std Deviation from MBAR
std_V_mbar = load('result_sigmaV.txt')*A32m3;
std_V2_mbar = load('result_sigmaV2.txt')*A32m3*A32m3;
std_H_mbar = load('result_sigmaH.txt')*kcalmol2J;
std_U_mbar = load('result_sigmaU.txt')*kcalmol2J;
std_H2_mbar = load('result_sigmaH2.txt')*kcalmol2J*kcalmol2J;
std_VH_mbar = load('result_sigmaVH.txt')*kcalmol2J*A32m3;

% Import NIST data
format long
filename = 'fluid.txt';
nist_data = load(filename);

%% Properties from NPT simulations
%%cp = (mean_H2-mean_H.^2)./(kb*mean_T.^2);
%cp = (mean_H2-mean_H.^2)./((n_atoms/n_a)*(kb*mean_T.^2)); %cp total
%Cpr = (mean_H2-mean_H.^2)./((n_atoms/n_a)*(kb*mean_T.^2));
Cpr = (mean_UHconf-mean_U.*mean_Hconf)./((n_atoms/n_a)*(kb*Ti.^2))+(mean_VHconf-mean_V.*mean_Hconf)./((n_atoms/n_a)*(kb*Ti.^2))+(19.5)*4.184; % (19.5-1.98720425864083)
Cprid = (mean_Hconf2-mean_Hconf.^2)./((n_atoms/n_a)*(kb*Ti.^2))+(19.5)*4.184; % cp conf + ideal
Cpt = (mean_H2-mean_H.^2)./((n_atoms/n_a)*(kb*mean_T.^2));% cp
Cp = Cprid;

aP = (mean_VH-mean_V.*mean_H)./(mean_V.*kb.*Ti.^2); %isobaric thermal expansion coeff
xT = (mean_V2 - mean_V.^2)./(mean_V.*kb.*Ti)*1e6; %isothermic compress %1e6 Pa to MPa
u_jt = (mean_V/(n_atoms/n_a)*1e6).*(Ti.*aP-1)./Cp; %1e6 Pa to MPa

Cv = Cp - Ti.*(mean_V/(n_atoms/n_a)*1e6).*aP.^2./xT; %1e6 Pa to MPa

rho = 1./(mean_V/(n_atoms/n_a)*1000); % mol/L

v_sound = sqrt(Cp./(Cv.*xT*1e-6.*rho*mole_mass)); % (m/s)

%% Properties from mbar
xT_mbar = (V2 - V.*V)./(kb*V.*T_mbar)*1e6;
aP_mbar = (VH - V.*H)./(kb*V.*(T_mbar.^2));
rho_mbar = 1./(V/(n_atoms/n_a)*1000);
Cp_mbar = (H2-H.*H)./((n_atoms/n_a)*kb*T_mbar.^2)+(19.5-1.98720425864083)*4.184;

%% Error bars
%std dev errors
%std_V2_mbar = sqrt(varV2);
%std_V_mbar = sqrt(varV);
% std_VH_mbar = sqrt(varVH);
% std_H_mbar = sqrt(varH);
dxT = 1./(kb.*Ti)*1e6 .* sqrt(std_V2_mbar.^2./mean_V.^2+std_V_mbar.^2.*((mean_V2+mean_V.^2)./mean_V.^2).^2);
daP = 1./(kb.*Ti.^2).*sqrt(std_VH_mbar.^2./mean_V.^2+std_H_mbar.^2+(mean_VH./mean_V.^2).^2.*std_V_mbar.^2);
dCp = 1./((n_atoms/n_a)*kb*Ti.^2).*sqrt(std_H2_mbar.^2+4*mean_H.^2.*std_H_mbar.^2);
dCv= sqrt(dCp.^2+((Ti.*aP.^2./xT*1e-6).^2).*std_V_mbar.^2+((2*Ti.*aP.*mean_V./xT*1e-6).^2.).*daP.^2+(Ti.*aP.^2.*mean_V./(xT*1e-6).^2).^2.*dxT.^2);
du_jt = 1e6*sqrt( ((Ti.*aP-1)./Cp).^2.*(std_V_mbar/(n_atoms/n_a)).^2 + (mean_V/(n_atoms/n_a).*Ti./Cp).^2.*daP.^2 + (mean_V/(n_atoms/n_a).*(Ti.*aP-1)./Cp.^2).^2.*dCp.^2 );
drho = rho.*std_V_mbar./mean_V;
dv_sound = sqrt(0.5*(Cp.*Cv.*xT*1e-6.*rho*mole_mass).^-1.*dCp.^2 + 0.5*Cp./(Cv.^3.*xT*1e-6.*rho*mole_mass).*dCv.^2 + 0.5*Cp./(Cv.*xT.^3*1e-6.*rho*mole_mass).*dxT.^2 + + 0.5*Cp./(Cv.*xT*1e-6.*rho.^3*mole_mass).*drho.^2);
%var errors
v = (mean_V2-mean_V.*mean_V).^2;
vh = (mean_VH-mean_V.*mean_H).^2;
uh = (mean_UH-mean_U.*mean_H).^2;

%dxT = xT.*sqrt((varV2) + (4.*mean_V.^2.*varV)./v) + (varV./mean_V.^2);
%daP = aP.*sqrt((varVH) + (mean_H.^2.*varV)+ (mean_V.^2.*varH)./vh) + (varV./mean_V.^2); 
drho = rho.*std_V_mbar./mean_V;

% P = 101325;
% a1 = (Cpr./(kb.*Ti.^2)-1000*kb)/kcalmol2J;
% dCpr = a1.*sqrt(((varUH) + (mean_H.^2.*varU) + (mean_U.^2.*varH))./uh) + P*(((varVH) + (mean_H.^2.*varV) + (mean_V.^2.*varH))./vh);
% dCp = dCpr*kcalmol2J+(19.5-1.98720425864083)*4.184; 

%CI - Confidencie Interval
% z = k/sqrt(N experiments)
z = 1.96/sqrt(1000000); %sugerir 5000 na razao e nao 1000

%% Plots
%Temperature convergence
% fig1 = figure(1);
% hold on
% cmean_T = zeros(20,1001);
% for i=1:20
%     cmean_T(i,:) = cummean(T(i,:));
%     plot(cmean_T(i,:))  
% end
% %ylim([270 360])
% hold off 

%% Ponto Critico
%Ponto Critico TraPPE/ Tc = 548.30; Pc = 47.67;
Tc_TraPPE = 548.30;
%Ponto Critico TraPPE/ Tc = 518.27; Pc = 46.75;
Tc_OPLS = 518.27;
%Ponto Critico TraPPE/ Tc = 562.27; Pc = 48.43;
Tc = 562.27;

%% PLOTS
fig2 =figure(2);
tt = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
xlabel (tt, 'Temperatura reduzida (T_r)','FontSize', 14, 'FontWeight', 'bold')

%% Xt
%%NIST DATA
xT_nist = [280, 8.35E-04;290,9.02E-04;300,9.77E-04;310,1.06E-03; 320,1.15E-03;330,1.26E-03;340,1.37E-03;350,1.50E-03]; % 360,1.02E+01 dado para gás
nexttile; 
plot (xT_nist(:,1)/Tc, xT_nist(:,2),'k-','LineWidth',3);
hold on
fit_xT_nist = fit(xT_nist(:,1),xT_nist(:,2),'poly2');
%plot(fit_xT_nist,'r-')
% MD DATA
errorbar(Ti/Tc_TraPPE,xT,z*dxT,'ok','MarkerFaceColor','k');
% plot(T_mbar,xT_mbar,'g-','LineWidth',1)
% MD DATA FIT
% fit_xT = fit(Ti,xT,'poly2');
% plot(fit_xT,'g-')
[p_xT,S_xT] = polyfit(Ti, xT, 2);
[xT_fit,xT_delta] = polyval(p_xT,Ti,S_xT);
plot(Ti/Tc_TraPPE,xT_fit, 'k-','LineWidth',1)
plot(Ti/Tc_TraPPE,xT_fit+2*xT_delta, 'k--',Ti/Tc_TraPPE,xT_fit-2*xT_delta, 'k--')
% PLOT DETAILS
%xlim([270 360])
ylim([.5e-3 1.9e-3])
ylabel ('X_{T} (MPa^{-1})','FontSize', 14, 'FontWeight', 'bold');
legend('NIST','TraPPE-EH G-ewald fixo','Curva de ajuste','Intervalo de Confiança (95%)','FontSize', 12')

%% a_P
% n = 370;
% Taux1 = linspace(nist_data(1,1), nist_data(n,1),n);
% rho_interp = polyfit(Taux1,nist_data(1:n,3),1);
% y = polyval(rho_interp,Taux1);
% dy = diff(y);
% dT = diff(Taux1);
% dydT = dy./dT;
% NIST_alpha_P1 = -dydT./nist_data(1:n-1,3).';
% %plot(Taux1(1:end-1),-dydT./nist_data(1:n-1,3).');
% %plot(Taux1,y, nist_data(:,1),nist_data(:,3));
% Taux2 = linspace(nist_data(n+1,1), nist_data(end,1),length(nist_data(:,1))-n);
% rho_interp = polyfit(Taux2,nist_data(n+1:end,3),1);
% y = polyval(rho_interp,Taux2);
% dy = diff(y);
% dT = diff(Taux2);
% dydT = dy./dT;
% NIST_alpha_P2 = -dydT./nist_data(n+1:end-1,3).';

n_trans = 370;  % transicao de fase ocorre na linha 370
j = 1;    % dy irá calcular com qual distância, derivada
alpha_P_NIST = zeros(n_trans-2*j,1);
T_NIST = zeros(n_trans-2*j,1);
for i=1+j:n_trans-j
    dy = nist_data(i+j,3)-nist_data(i-j,3);
    dt = nist_data(i+j,1)-nist_data(i-j,1);
    alpha_P_NIST(i-1) = -1./nist_data(i,3).*(dy./dt);
    T_NIST(i-1) = nist_data(i,1);
end
fit_aP_nist = fit(T_NIST,alpha_P_NIST,'poly2');

%fig3 = figure(3);
nexttile; 
%old plot
%plot (T_NIST,alpha_P_NIST,'r-', Taux2(1:end-1), NIST_alpha_P2,'r-',Ti,aP,'b-o','LineWidth',2);
% NIST DATA
plot (T_NIST/Tc,alpha_P_NIST,'k-','LineWidth',3);
hold on
% MD DATA
errorbar(Ti/Tc_TraPPE,aP,z*daP,'ok','MarkerFaceColor','k');
% MD DATA FIT
[p_aP,S_aP] = polyfit(Ti, aP, 2);
[aP_fit,aP_delta] = polyval(p_aP,Ti,S_aP);
plot(Ti/Tc_TraPPE,aP_fit, 'k-','LineWidth',1)
plot(Ti/Tc_TraPPE,aP_fit+2*aP_delta, 'k--',Ti/Tc_TraPPE,aP_fit-2*aP_delta, 'k--')
%fit_aP = fit(Ti,aP,'poly2');
%plot(fit_aP,'g-')
%plot (T_mbar,aP_mbar,'g-','LineWidth',1);
% PLOT DETAILS
%xlabel ('T_r', 'FontSize', 14, 'FontWeight', 'bold');
%xlim([270 360])
ylim([1e-3 1.9e-3])
ylabel ('a_{P} (K^{-1})', 'FontSize', 14, 'FontWeight', 'bold');
% legend ('NIST','TraPPE','FIT POLY2','95% Predction Interval','FontSize', 12)

%% Density
% rho_c = 3.9;
% rho_c_TraPPE = 3.8;
% %fig4 = figure(4);
% % NIST DATA
% nexttile; plot(nist_data(:,1)/Tc, nist_data(:,3)/rho_c,'k-','LineWidth',3)
% hold on
% % MD DATA
% errorbar(Ti/Tc_TraPPE,rho/rho_c_TraPPE ,z*drho,'ok','MarkerFaceColor','k');
% % plot(T_mbar, rho_mbar, 'g-')
% % MD DATA FIT
% [p_rho,S_rho] = polyfit(Ti, rho, 1);
% [rho_fit,rho_delta] = polyval(p_rho,Ti,S_rho);
% plot(Ti/Tc_TraPPE,rho_fit/rho_c_TraPPE , 'k-','LineWidth',1)
% %plot(Ti,rho_fit+2*aP_delta, 'g--',Ti,rho_fit-2*aP_delta, 'g--')
% % PLOT DETAILS
% %xlim([270 390])
% % ylim([10 12])
% %legend ('NIST','TraPPE-EH Greedy Search', 'TraPPE-EH','FontSize', 14)
% %plot(data(:,1), data(:,3))
% legend ('NIST','TraPPE','LINEAR FIT','FontSize', 12)
% ylabel ('Density (mol/L)','FontSize', 14, 'FontWeight', 'bold');
% %xlabel ('T_r', 'FontSize', 14, 'FontWeight', 'bold');
% hold off

%% C_P
%fig5 = figure(5);
nexttile; 
% NIST DATA
plot(nist_data(:,1)/Tc, nist_data(:,9),'k-','LineWidth',3);
hold on
% MD DATA
errorbar(Ti/Tc_TraPPE,Cp,z*dCp,'ok','MarkerFaceColor','k');
% MD DATA FIT
[p_Cp,S_Cp] = polyfit(Ti, Cp, 2);
[Cp_fit,Cp_delta] = polyval(p_Cp,Ti,S_Cp);
plot(Ti/Tc_TraPPE,Cp_fit, 'k-','LineWidth',1)
plot(Ti/Tc_TraPPE,Cp_fit+2*Cp_delta, 'k--',Ti/Tc_TraPPE,Cp_fit-2*Cp_delta, 'k--')
%fit_Cp = fit(Ti,Cp,'poly2');
fit_Cp_nist = fit(nist_data(1:n_trans,1),nist_data(1:n_trans,9),'poly2');
% PLOT DETAILS
% plot(fit_Cp_nist,'g-')
% plot(T_mbar,Cp_mbar,'g-','LineWidth',1)
% plot(fit_Cp, 'g-')
% legend ('NIST','TraPPE','FIT POLY2','95% Predction Interval','FontSize', 12)
ylabel('C_{P} (J/mol.K)','FontSize', 14, 'FontWeight', 'bold');
%xlim([270 360])
% ylim([125 160])

%% C_V
nexttile; 
% NIST DATA
plot(nist_data(:,1)/Tc, nist_data(:,8),'k-','LineWidth',3);
hold on
% MD DATA
errorbar(Ti/Tc_TraPPE,Cv,z*dCv,'ok','MarkerFaceColor','k');
% MD DATA FIT
[p_Cv,S_Cv] = polyfit(Ti, Cv, 2);
[Cv_fit,Cv_delta] = polyval(p_Cv,Ti,S_Cv);
plot(Ti/Tc_TraPPE,Cv_fit, 'k-','LineWidth',1)
plot(Ti/Tc_TraPPE,Cv_fit+2*Cv_delta, 'k--',Ti/Tc_TraPPE,Cv_fit-2*Cv_delta, 'k--')
%fit_Cv = fit(Ti,Cv,'poly2');
fit_Cv_nist = fit(nist_data(1:n_trans,1),nist_data(1:n_trans,8),'poly2');
% PLOT DETAILS
% legend ('NIST','TraPPE','FIT POLY2','95% Predction Interval','FontSize', 12)
ylabel('C_{V} (J/mol.K)','FontSize', 14, 'FontWeight', 'bold');
%xlim([270 360])
ylim([75 110])

%% joule thompson
nexttile;
% NIST DATA
u_jt_NIST = nist_data(1:n_trans,11)*5/9/0.101325;
plot(nist_data(1:n_trans,1)/Tc, u_jt_NIST, 'k-','LineWidth',3);
% MD DATA
hold on
errorbar(Ti/Tc_TraPPE,u_jt,z*du_jt,'ok','MarkerFaceColor','k');
% MD DATA FIT
[p_u_jt,S_u_jt] = polyfit(Ti, u_jt, 2);
[u_jt_fit,u_jt_delta] = polyval(p_u_jt,Ti,S_u_jt);
plot(Ti/Tc_TraPPE,u_jt_fit, 'k-','LineWidth',1)
plot(Ti/Tc_TraPPE,u_jt_fit+2*u_jt_delta, 'k--',Ti/Tc_TraPPE,u_jt_fit-2*u_jt_delta, 'k--')
%fit_jt = fit(Ti,u_jt,'poly2');
%plot(fit_jt,'g-')
% PLOT DETAILS
% legend ('NIST','TraPPE','FIT POLY2','95% Predction Interval','FontSize', 12)
ylabel('u_{JT} (K/MPa)','FontSize', 14, 'FontWeight', 'bold');
ylim([-.5 -0.2])
% xlim([270 360])

%% Sound Speed
nexttile;
% NIST DATA
v_sound_NIST = (nist_data(1:n_trans,10));
plot(nist_data(1:n_trans,1)/Tc, v_sound_NIST, 'k-', 'LineWidth',3);
[p_v_sound_nist,S_v_sound_nist] = polyfit(nist_data(1:n_trans,1), v_sound_NIST, 1);
hold on
% MD DATA
errorbar(Ti/Tc_TraPPE,v_sound,z*dv_sound,'ok','MarkerFaceColor','k');
% MD DATA FIT
[p_v_sound,S_v_sound] = polyfit(Ti, v_sound, 1);
[v_sound_fit,v_sound_delta] = polyval(p_v_sound,Ti,S_v_sound);
plot(Ti/Tc_TraPPE,v_sound_fit, 'k-','LineWidth',1)
plot(Ti/Tc_TraPPE,v_sound_fit+2*v_sound_delta, 'k--',Ti/Tc_TraPPE,v_sound_fit-2*v_sound_delta, 'k--')
% fit_v_sound = fit(Ti,v_sound, 'poly1');
% plot(fit_v_sound,'g-');
% legend ('NIST','TraPPE','FIT POLY1','FontSize', 12)
ylabel('v_{som}(m/s)','FontSize', 14, 'FontWeight', 'bold');

% %% Enthalpy 
% nexttile;
Enthalpy = mean_H*(n_a/n_atoms)/1000; %kJ/mol 
% plot (nist_data(1:n_trans,1)/Tc, nist_data(1:n_trans,6),'r-',Ti/Tc_TraPPE,Enthalpy,'bo','LineWidth',2);
[p_H_nist,S_H_nist] = polyfit(nist_data(1:n_trans,1), nist_data(1:n_trans,6), 1);
[p_H,S_H] = polyfit(Ti, Enthalpy, 1);
[H_fit,H_delta] = polyval(p_H,Ti,S_H);
% hold on
% % errorbar(Ti,Enthalpy,z*std_H_mbar*(n_a/n_atoms)/1000,'ob','MarkerFaceColor','b');
% % xlim([270 390])
% ylim([-40 40])
% xlabel ('T','FontSize', 14, 'FontWeight', 'bold')
% ylabel ('Enthalpy (kJ/mol)', 'FontSize', 14, 'FontWeight', 'bold');

%% Volume
% nexttile;
% Volume = 1./rho; %kJ/mol 
% plot (Ti,Volume,'bo','LineWidth',2);
% hold on
% fit_V = fit(Ti,Volume,'poly3');
% V_coeffs = coeffvalues(fit_V);
% plot(fit_V,'g-')
% xlabel ('T','FontSize', 14, 'FontWeight', 'bold')
% ylabel ('Enthalpy (kJ/mol)', 'FontSize', 14, 'FontWeight', 'bold');

%% Internal Energy
%fig6 = figure(6);
% nexttile;
% IntEnergy = mean_Etotal*(n_a/n_atoms)/1000; %kJ/mol
% plot (nist_data(1:n_trans,1), nist_data(1:n_trans,6),'r-', Ti,IntEnergy,'bo','LineWidth',2);
% hold on
% errorbar(Ti,IntEnergy,z*std_U_mbar*(n_a/n_atoms)/1000,'ob','MarkerFaceColor','b');
% xlim([270 390])
% ylim([-30 40])
% xlabel ('T','FontSize', 14, 'FontWeight', 'bold')
% ylabel ('Internal Energy (kJ/mol)', 'FontSize', 14, 'FontWeight', 'bold');

%% Error PLOTS
Tplot = linspace(280,350);
%xT
coeff_xT_nist = coeffvalues(fit_xT_nist);
coeff_xT = p_xT;
ERROR_xT = coeff_xT_nist(1).*Tplot.^2+coeff_xT_nist(2).*Tplot+coeff_xT_nist(3)-coeff_xT(1).*Tplot.^2-coeff_xT(2).*Tplot-coeff_xT(3);
ABS_xT = 100*abs(ERROR_xT./(coeff_xT_nist(1).*Tplot.^2+coeff_xT_nist(2).*Tplot+coeff_xT_nist(3)));
%aP
coeff_aP_nist = coeffvalues(fit_aP_nist);
coeff_aP = p_aP;
ERROR_aP = abs(coeff_aP_nist(1).*Tplot.^2+coeff_aP_nist(2).*Tplot+coeff_aP_nist(3)-coeff_aP(1).*Tplot.^2-coeff_aP(2).*Tplot-coeff_aP(3));
ABS_aP = 100*abs(ERROR_aP./(coeff_aP_nist(1).*Tplot.^2+coeff_aP_nist(2).*Tplot+coeff_aP_nist(3)));
%densidade
fit_rho = fit(Ti,rho,'poly1');
fit_rho_nist = fit(nist_data(1:n_trans,1), nist_data(1:n_trans,3),'poly1');
ERROR_rho = coeffvalues(fit_rho_nist)-coeffvalues(fit_rho);
ERROR_rho = ERROR_rho(1).*Tplot+ERROR_rho(2);
fit_rho_nist = coeffvalues(fit_rho_nist);
ABS_rho = 100*abs(ERROR_rho./(fit_rho_nist(1).*Tplot+fit_rho_nist(2)));
%aad_rho = mad(ABS_rho,0,'all');
%u_jt
fit_jt_nist = fit(nist_data(1:n_trans,1), u_jt_NIST,'poly2');
coeff_jt_nist = coeffvalues(fit_jt_nist);
coeff_jt = p_u_jt;
ERROR_jt = coeff_jt_nist(1).*Tplot.^2+coeff_jt_nist(2).*Tplot+coeff_jt_nist(3)-coeff_jt(1).*Tplot.^2-coeff_jt(2).*Tplot-coeff_jt(3);
ABS_jt = 100*abs(ERROR_jt./(coeff_jt_nist(1).*Tplot.^2+coeff_jt_nist(2).*Tplot+coeff_jt_nist(3)));
%Cp
coeff_cp_nist = coeffvalues(fit_Cp_nist);
coeff_cp = p_Cp;
ERROR_cp = coeff_cp_nist(1).*Tplot.^2+coeff_cp_nist(2).*Tplot+coeff_cp_nist(3)-coeff_cp(1).*Tplot.^2-coeff_cp(2).*Tplot-coeff_cp(3);
ABS_cp = 100*abs(ERROR_cp./(coeff_cp_nist(1).*Tplot.^2+coeff_cp_nist(2).*Tplot+coeff_cp_nist(3)));
%Cv
coeff_Cv_nist = coeffvalues(fit_Cv_nist);
coeff_Cv = p_Cv;
ERROR_Cv = coeff_Cv_nist(1).*Tplot.^2+coeff_Cv_nist(2).*Tplot+coeff_Cv_nist(3)-coeff_Cv(1).*Tplot.^2-coeff_Cv(2).*Tplot-coeff_Cv(3);
ABS_Cv = 100*abs(ERROR_Cv./(coeff_Cv_nist(1).*Tplot.^2+coeff_Cv_nist(2).*Tplot+coeff_Cv_nist(3)));
% Sound Speed
coeff_v_sound = p_v_sound;
ERROR_v_sound = p_v_sound_nist(1).*Tplot + p_v_sound_nist(2) - (p_v_sound(1).*Tplot+p_v_sound(2));
ABS_v_sound = 100*abs(ERROR_v_sound./(p_v_sound_nist(1).*Tplot + p_v_sound_nist(2)));
% Enthalpy
coeff_H= p_H;
ERROR_H = p_H_nist(1).*Tplot + p_H_nist(2) - (p_H(1).*Tplot+p_H(2));
ABS_H = 100*abs(ERROR_H./(p_H_nist(1).*Tplot + p_H_nist(2)));

%% AAD
AAD = [sum(ABS_rho,'all')/length(ABS_rho),
        sum(ABS_xT,'all')/length(ABS_xT),
        sum(ABS_aP,'all')/length(ABS_aP),
        sum(ABS_jt,'all')/length(ABS_jt),
        sum(ABS_v_sound,'all')/length(ABS_v_sound),
        sum(ABS_cp,'all')/length(ABS_cp),
        sum(ABS_Cv,'all')/length(ABS_Cv),
        sum(ABS_H,'all')/length(ABS_H)];
%% ERRORS plot
fig4 = figure(4);
tiledlayout(1,2);
% nexttile;
% PLOT ERROR NORMAL SCALE
% plot(Tplot,ABS_rho,'k-',Tplot,ABS_xT,'k--', Tplot,ABS_aP,'k-.', Tplot,ABS_jt,'k:','LineWidth',3)
% hold on
% plot(Tplot,ABS_cp,'k-o',Tplot,ABS_Cv, 'k-+','LineWidth',2)
% plot(Tplot,ABS_v_sound, 'k-')

nexttile;
semilogy(Tplot,ABS_cp,'ko',Tplot,ABS_Cv, 'k+','LineWidth',2)
hold on
semilogy(Tplot,ABS_v_sound, 'k-.','LineWidth',2)
ylim([1e-2 1e2])
ylabel ('Relative Error (%)', 'FontSize', 14, 'FontWeight', 'bold');
legend('C_P','C_V','v_{sound}')
grid on

nexttile;
semilogy(Tplot,ABS_rho,'k-',Tplot,ABS_xT,'k--', Tplot,ABS_aP,'k-.', Tplot,ABS_jt,'k:','LineWidth',3)
hold on
grid on
ylim([1e-2 1e2])
legend('rho','x_T','a_P','u_{JT}')
ylabel ('Relative Error (%)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel ('Temperature','FontSize', 14, 'FontWeight', 'bold')