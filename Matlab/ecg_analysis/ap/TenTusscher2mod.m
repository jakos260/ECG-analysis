function [t,V]=TenTusscher2mod(HT, STOPTIME, Stim_I, Stim_T, ISO, Stim_Int, CT, phases_mod)
%% TENTUSSCHER2
%     Ten Tusscher model, 2nd version, described in the article: KHWJ ten
%     Tusscher and AV Panfilov, "Alternans and spiral breakup in a human
%     ventricular tissue model", Am J Physiol Heart Circ Physiol, 291:
%     H1088-H1100, 2006. (Note that the parameter values for Vc, Vsr, Vss,
%     Vrel, and k4 as in the article, were later corrected by Ten Tusscher)
%
% VERSION
%     20100922V0.2   (22 Sept. 2010)
%     This version is a MatLab rewrite of the C++ version downloaded from
%     Ten Tusscher's website: http://www-binf.bio.uu.nl/khwjtuss/
%
% INPUT VARIABLES
%     HT       : Timestep (ms)
%     STOPTIME : Duration of the simulation (ms)
%     Stim_I   : Array with the stimulus strengths of monophasic square
%                stimulus pulses (A/F). One entry per pulse.
%     Stim_T   : Array with the respective pulse durations (ms)
%     ISO      : Degree of sympathetic innervation. Float between 0 and 1
%     Stim_Int : Array with interpulse intervals (ms). The interpulse
%                interval is the interval after the start of a pulse to the
%                start of the next pulse or to the end of the simulation.
%     CT       : Cell type: 1=epicardial, 2=mid-myocardial, 3=endocardial
%     phase_mod   : Array of 3 factors modifing 3 phases currents
%
% OUTPUT VARIABLES
%     t : Array of M points in time for each sample.
%     V : MxN array containing the output variables. In the current version
%         N=1 and V contains only the action potential at M time points as
%         defined in t.
%
% AUTHOR (translator)
%     Johannes Jan Struijk
%     jjs@hst.aau.dk
%     Laboratory for Cardiotechnology
%     Dept. Health Science and Technology
%     Aalborg University
%     Denmark
%

%% Flag to choose type of cell
switch CT
    case 1, CellType = 'Epicardial';
    case 2, CellType = 'Mid-myocardial';
    case 3, CellType = 'Endocardial';
end

%% Parameters of beta-adrenergic response

%% from Gong Article 
beta_GCal = 1+(1.9-1)*ISO;
beta_f_shift_ICal = 8*ISO;
beta_d_shift_ICal = 12*ISO;
beta_GKs = 1+(8-1)*ISO;
beta_GNa = 1+(2.7-1)*ISO;
beta_h_INa_shift = 5*ISO;
beta_kNai = 1+(0.7-1)*ISO;
Irel_beta = 1+(2.5-1)*ISO;
Iup_scale = 1+(0.54-1)*ISO;
IKs_beta_shift = 10*ISO;
% shifting parameters for IKs is depending on 

%Change from Zhang et al. 2026
IpK_beta = 1+(2-1)*ISO;

%% Electrophysiological parameters

Ko  =   5.4;      % Extracellular K concentration (mM)
Cao =   2.0;      % Extracellualr Ca concentration (mM)
Nao = 140.0;      % Extracellular Na concentration (mM)

Vc  = 0.016404;   % Volume of the cytoplasm (mm^3) 
Vsr = 0.001094;   % Volum of sarcoplasmatic reticulum (mm^3)
Vss = 0.00005468; % Volume of the subspace (mm^3)

Capacitance = 0.185; % Cellular capacitance (nF??)

Bufc   = 0.2;     % Cytoplasmic Ca buffer concentration (mM)
Kbufc  = 0.001;   % Cai, half-saturation constant for cytopl. buffer (mM)
Bufsr  = 10.0;    % Sarcoplasmic Ca buffer concentration (mM)
Kbufsr = 0.3;     % CaSR half-saturation constant for SR buffer (mM)
Bufss  = 0.4;     % Subspace Ca buffer concentration (mM)
Kbufss = 0.00025; % CaSS half-saturation constant for SS buffer (mM)

Vmaxup = 0.000425;% Maximal Iup conductance (mM/ms)
Kup    = 0.00025; % Half-saturation constant of Iup (mM) 
Vrel   = 0.102;   % Maximal Irel conductance (mM/ms)
k1_    = 0.15;    % Irel transition rate R to O and RI to I (1/mM^2*ms)
k2_    = 0.045;   % Irel transition rate O to I and R to RI (1/mM*ms)
k3     = 0.06;    % Irel transition rate O to R and I to RI (1/ms)
k4     = 0.005;   % Irel transition rate I to O and RI to I (1/ms)
EC     = 1.5;     % CaSR half-saturation constant of kCaSR (mM)
maxsr  = 2.5;     % Maximum value of kCaSR (dimensionless)
minsr  = 1.0;     % Minimum value of kCaSR (dimensionless)
Vleak  = 0.00036; % Maximal Ileak conductance (mM/ms)
Vxfer  = 0.0038;  % Maximal Ixfer conductance (mM/ms)

R = 8314.472;     % Gas constant (mJ/K*mol)
F = 96485.3415;   % Faraday constant (C/mol)
T = 310;          % Temperature (K)
RTonF = R*T/F;


%% Parameters for ion currents

GKr=0.046;

switch CellType
    case 'Epicardial',     GKr = GKr*1.3; % Maximal IKs conductance (nS/pF)
    case 'Mid-myocardial', GKr = GKr*0.8; % "
    case 'Endocardial',    GKr = 0.046; % "
end

GK1=0.1908;
%GK1 = 5.405;        % Maximal IK1 conductance (nS/pF)
% GK1_0 = 0.6821;
% GK1_old = 5.405;
switch CellType
    case 'Epicardial',     GK1=GK1*1.2; % Maximal IKs conductance (nS/pF)
    case 'Mid-myocardial', GK1=GK1*1.3; % "
    case 'Endocardial',    GK1=0.1908;  % "
end

switch CellType
    case 'Epicardial',     GKs = 0.392*beta_GKs; % Maximal IKs conductance (nS/pF)
    case 'Mid-myocardial', GKs = 0.098*beta_GKs; % "
    case 'Endocardial',    GKs = 0.392*beta_GKs; % "
end
%GK1 = 5.405;        % Maximal IK1 conductance (nS/pF) 
switch CellType
    case 'Epicardial',     Gto = 0.294; % Maximal Ito conductance (nS/pF)
    case 'Mid-myocardial', Gto = 0.294; % "
    case 'Endocardial',    Gto = 0.073; % "
end
%GNa   = 14.838*beta_GNa;     % Maximal INa conductance (nS/pF)
GNa   = 75*beta_GNa; %14.838;     % Maximal INa conductance (nS/pF)
switch CellType
    case 'Epicardial',     GNaL = 0.0075*0.6; % Maximal INaL conductance (nS/pF)
    case 'Mid-myocardial', GNaL = 0.0075; % "
    case 'Endocardial',    GNaL = 0.0075; % "
end

GbNa  = 0.00029;    % Maximal background-INa conductance (nS/pF)
KmK   = 1.0;
KmNa  = 40.0*beta_kNai;
knak  = 2.724;
pKNa  = 0.03;
GCaL  = 0.0000398*beta_GCal;
GbCa  = 0.000592;
knaca = 1000;
KmNai = 87.5;
KmCa  = 1.38;
ksat  = 0.1;
n     = 0.35;       % Voltage dependence parameter of INaCa
GpCa  = 0.1238;
KpCa  = 0.0005;
GpK   = 0.0146;

% % magnesium parameters
% p_GK1 = 0.8838;
% s_Mg = 1.0648;
% SPMc = 0.0014613;

% % magnesium parameters (constants)
% bufmg  = 5.67;   % mM
% KMgbuf = 0.174;  % mM 
% Mgtot  = 5.0;    % mM
% 
% % compute equilibrium free Mg2+ from quadratic
% b_Mgi = bufmg - Mgtot + KMgbuf;
% c_Mgi = KMgbuf*Mgtot;
% Mgi   = (-b_Mgi + sqrt(b_Mgi^2 + 4*c_Mgi)) / 2;   % ≈ 0.6 mM


%% Initial conditions  for state variables
inverseVcF   = 1/(Vc*F);
inverseVcF2  = inverseVcF/2;
inverseVssF2 = 1/(2*Vss*F);
svolt  = -86.2;
Cai    = 0.00007;
CaSR   = 1.3;
CaSS   = 0.0007;
Nai    = 7.67;
Ki     = 138.3;
%sh     = 0.75;
%sj     = 0.75;
rkr   = 0.0;
xr   = 1.0;
sxs    = 0.0;
sr     = 0.0;
ss     = 1.0;
sd     = 0.0;
sf     = 1.0;
sf2    = 1.0;
sfcass = 1.0;
sRR    = 1.0;
sOO    = 0.0;
xrs    = 0.5;
xrf    = 0.0;
rk1    = 0.5;
xk1    = 0.5;
A_H_f  = 0.99;
A_H_s  = 1.0-A_H_f;
shf    = 0.75;
shs    = 0.75;
sj     = 0.75;
smL    = 0;
shL    = 1;
sm     = 0.0;

% End of Electrophysiological parameters

%% Integration parameters

% HT = 0.02;                       % timestep (ms)
% STOPTIME = 1000 + sum(Stim_Int); % Duration of the simulation (ms)

% End of Integration parameters

LEN = ceil(STOPTIME/HT)+1;
Vm = zeros(LEN,1);         % Membrane potential (mV)
SINaCa = zeros(LEN,1);     % INaCa to be saved
SIKr = zeros(LEN,1);
SIKs = zeros(LEN,1);
SIK1 = zeros(LEN,1);
SIto = zeros(LEN,1);
SINa = zeros(LEN,1);
SINaL = zeros(LEN,1);
SIbNa = zeros(LEN,1);
SICaL = zeros(LEN,1);
SIbCa = zeros(LEN,1);
SINaK = zeros(LEN,1);
SIpCa = zeros(LEN,1);
SIpK = zeros(LEN,1);
SKi = zeros(LEN,1);
SCai = zeros(LEN,1);
SNai = zeros(LEN,1);
sItot = zeros(LEN,1);
TimeArray = zeros(LEN,1);  % Timestamp for samples (ms)
savetime = 0;              % Time to save parameters (every appr. 1ms)
savecount = 1;             % Nr of times parameters have been saved

% End of initial values

%% THE MAIN LOOP STARTS HERE

Istim = zeros(LEN,1);
pulseStart = 1;
for pulseNr = 1:length(Stim_I)
    pulseStart = pulseStart + ceil(Stim_Int(pulseNr)/HT);
    pulseStop  = pulseStart + floor(Stim_T(pulseNr)/HT);
    Istim(pulseStart:pulseStop) = ones(pulseStop-pulseStart+1,1)*Stim_I(pulseNr);
end

time = 0;
while time < STOPTIME

    time = time + HT;

    % Parameters needed to calculate ion currents
    EK  = RTonF * log(Ko/Ki);
    ENa = RTonF * log(Nao/Nai);
    EKs = RTonF * log((Ko+pKNa*Nao)/(Ki+pKNa*Nai));
    ECa = 0.5*RTonF*log(Cao/Cai);
    %AK1 = 0.1/(1+exp(0.06*(svolt-EK-200)));
    %BK1 = (3*exp(0.0002*(svolt-EK+100)) + exp(0.1*(svolt-EK-10))) / ...
    %       (1+exp(-0.5*(svolt-EK)));
    %rec_iK1  = AK1/(AK1+BK1);
    rec_iNaK = 1 / ...
                (1+0.1245*exp(-0.1*svolt/RTonF)+0.0353*exp(-svolt/RTonF));
    rec_ipK = 1 / (1+exp((25-svolt)/5.98));
    sh=A_H_f*shf+A_H_s*shs;


    % % use Mgi as free Mg2+ (constant or slowly varying if you later add dynamics)
    % KI_Mg   = 2.8*exp(-(svolt - s_Mg*EKs)/180);
    % KB_Mg   = 0.45*exp(-(svolt - s_Mg*EKs)/20);
    % Kd1_spm = 0.0007*exp(-(svolt - s_Mg*EKs + 8*Mgi)/4 - 8);
    % Kd2_spm = 0.04*exp(-(svolt - s_Mg*EKs)/9.1);
    % psi_Mg  = 1 + Mgi/KB_Mg;
    % r1_inf  = psi_Mg^2 / ((SPMc/Kd1_spm) + (Mgi/KI_Mg) + psi_Mg^3);
    % r2_inf  = 1 / (1 + (SPMc/Kd2_spm));
    % G_Gmax_ratio = p_GK1*r1_inf + (1 - p_GK1)*r2_inf;
    % 
    % a_G = 1/35;
    % b_G = (-55)/7;
    %GK1   = GK1_0*(a_G*T + b_G)*sqrt(Ko/5.4);
    
    
    % Calculate ion currents
    INa  = (GNa*(sm^3)*sh*sj*(svolt-ENa));
    INaL = GNaL*(svolt-ENa)*smL*shL;
    svolt_etc = 2*(svolt-15)/RTonF;
    ICaL = (GCaL*sd*sf*sf2*sfcass*2*F*(svolt_etc) * ...
           (0.25*exp(svolt_etc)*CaSS-Cao) / (exp(svolt_etc)-1));
    Ito = (Gto*sr*ss*(svolt-EK));
    IKr = (GKr*sqrt(Ko/5.4)*xr*rkr*(svolt-EK));
    IKs = (GKs*sxs*sxs*(svolt-EKs));
    %IK1 = GK1*rec_iK1*(svolt-EK);
    IK1= (GK1*sqrt(Ko)*rk1*xk1*((svolt)-EK)); %Ohara Rudy
    %IK1 = GK1*G_Gmax_ratio*(svolt-EK); % FINK
    INaCa = knaca*(1/((KmNai^3)+(Nao^3)))*(1/(KmCa+Cao))* ...
            (1/(1+ksat*exp((n-1)* svolt/RTonF)))* ...
            (exp(n*svolt/RTonF)*(Nai^3)*Cao-exp((n-1)*svolt/RTonF)* ...
            (Nao^3)*Cai*2.5);
    INaK = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
%   INaK = 0.5*knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa = GpCa*Cai/(KpCa+Cai);
    IpK  = ((GpK*rec_ipK*(svolt-EK))*IpK_beta);
    IbNa = GbNa*(svolt-ENa);
    IbCa = GbCa*(svolt-ECa);

    % Phases modifications
    phase1_I_mod = phases_mod(1);
    phase2_I_mod = phases_mod(2);
    phase3_I_mod = phases_mod(3);

    Ito     = Ito   * phase1_I_mod;
    ICaL    = ICaL  * phase2_I_mod;
    IKr     = IKr   * phase3_I_mod;
    IKs     = IKs   * phase3_I_mod;
    
    % Total current is sum of all the above
    sItot = IKr+IKs+IK1+Ito+INa+INaL+IbNa+ICaL+IbCa+INaK+INaCa+IpCa+IpK+Istim(savecount);
    
    % Update concentrations
    kCaSR = maxsr-(maxsr-minsr)/(1+(EC/CaSR)^2);
    k1    = k1_/kCaSR;
    k2    = k2_*kCaSR;
    dRR   = k4*(1-sRR)-k2*CaSS*sRR;
    sRR   = sRR+HT*dRR;
    sOO   = k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);
    
    Irel  = (Vrel*sOO*(CaSR-CaSS))*Irel_beta;
    Ileak = Vleak*(CaSR-Cai);
    Iup   = Vmaxup/((1+(Kup/Cai)^2)*Iup_scale);
    Ixfer = Vxfer*(CaSS-Cai);
    
    CaCSQN = Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR  = HT*(Iup-Irel-Ileak);
    bjsr   = Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr   = Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR   = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
    
    CaSSBuf = Bufss*CaSS/(CaSS+Kbufss);
    dCaSS   = HT*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inverseVssF2*...
               Capacitance));
    bcss    = Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
    ccss    = Kbufss*(CaSSBuf+dCaSS+CaSS);
    CaSS    = (sqrt(bcss*bcss+4*ccss)-bcss)/2;

    if CaSS > 0.025
        CaSS = 0.025;
    else
        CaSS = CaSS;
    end
    
    CaBuf = Bufc*Cai/(Kbufc+Cai);
    dCai  = HT*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*Capacitance)- ...
             (Iup-Ileak)*(Vsr/Vc)+Ixfer);
    bc    = Bufc-CaBuf-dCai-Cai+Kbufc;
    cc    = Kbufc*(CaBuf+dCai+Cai);
    Cai   = (sqrt(bc*bc+4*cc)-bc)/2;
    
    dNai = -(INa+INaL+IbNa+3*INaK+3*INaCa)*inverseVcF*Capacitance;
    Nai  = Nai + HT*dNai;
    
    dKi = -(Istim(savecount)+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*Capacitance;
    Ki  = Ki + HT*dKi;
        
    % Calculate steady state values and time constants (gating parameters)
    % AM    = 1/(1+exp((-60-svolt)/5));
    % BM    = 0.1/(1+exp((svolt+35)/5))+0.1/(1+exp((svolt-50)/200));
    % Tau_M = AM*BM;
    % M_Inf = 1/((1+exp((-56.86-svolt)/9.03))^2);
    % if svolt>-40
    %     AH_1  = 0;
    %     BH_1  = 0.77/(0.13*(1+exp(-(svolt+10.66)/11.1)));
    %     Tau_H =1/(AH_1+BH_1);
    % else
    %     AH_2  = (0.057*exp(-(svolt+80)/6.8));
    %     BH_2  = 2.7*exp(0.079*svolt)+3.1e5*exp(0.3485*svolt);
    %     Tau_H = 1/(AH_2+BH_2);
    % end
    % H_Inf = 1/((1+exp(svolt+71.55+beta_h_INa_shift)/7.43))^2;
    % if svolt>-40
    %     AJ_1  = 0;
    %     BJ_1  = (0.6*exp((0.057)*svolt)/(1+exp(-0.1*(svolt+32))));
    %     Tau_J = 1/(AJ_1+BJ_1);
    % else
    %     AJ_2  = (((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)* ...
    %              exp(-0.04391*svolt))*(svolt+37.78)/ ...
    %              (1+exp(0.311*(svolt+79.23))));
    %     BJ_2  = (0.02424*exp(-0.01052*svolt)/ ...
    %              (1+exp(-0.1378*(svolt+40.14))));
    %     Tau_J = 1/(AJ_2+BJ_2);
    % end
    % J_Inf = 1/((1+exp(svolt+71.55)/7.43))^2;

    M_Inf   = 1.0/(1.0+exp((-(svolt+39.57))/9.871));
    Tau_M   = 1.0/(6.765*exp((svolt+11.64)/34.77)+8.552*exp(-(svolt+77.42)/5.955));
    H_Inf   = 1.0/(1+exp((svolt+82.90+beta_h_INa_shift)/6.086));
    Tau_H_f = 1.0/(1.432e-5*exp(-(svolt+1.196)/6.285)+6.149*exp((svolt+0.5096)/20.27));
    Tau_H_s = 1.0/(0.009794*exp(-(svolt+17.95)/28.05)+0.3343*exp((svolt+5.730)/56.66));
    J_Inf   = H_Inf;
    Tau_J   = 2.038+1.0/(0.02136*exp(-(svolt+100.6)/8.281)+0.3052*exp((svolt+0.9941)/38.45));
    
    ML_Inf  = 1.0/(1.0+exp((-(svolt+42.85))/5.264));
    Tau_ML  = Tau_M;
    HL_Inf  = 1.0/(1.0+exp((svolt+87.61)/7.488));
    Tau_HL  = 200.0;
    
    xrss=1.0/(1.0+exp((-(svolt+8.337))/6.789));
    txrf=12.98+1.0/(0.3652*exp((svolt-31.66)/3.869)+4.123e-5*exp((-(svolt-47.78))/20.38));
    txrs=1.865+1.0/(0.06629*exp((svolt-34.70)/7.355)+1.128e-5*exp((-(svolt-29.74))/25.94));
    Axrf=1.0/(1.0+exp((svolt+54.81)/38.21));
    Axrs=1.0-Axrf;
    xr=Axrf*xrf+Axrs*xrs;
    rkr=1.0/(1.0+exp((svolt+55.0)/75.0))*1.0/(1.0+exp((svolt-10.0)/30.0));

    xk1ss=1.0/(1.0+exp(-(svolt+2.5538*Ko+144.59)/(1.5692*Ko+3.8115)));
    txk1=122.2/(exp((-(svolt+127.2))/20.36)+exp((svolt+236.8)/69.33));
    rk1=1.0/(1.0+exp((svolt+105.8-2.6*Ko)/9.493));
    
    Xs_Inf  = 1/(1+exp((-5-svolt)/14));
    Axs     = 1400/(sqrt(1+exp((5-svolt)/6)));
    Axs_beta     = 1400/(sqrt(1+exp((5-svolt+IKs_beta_shift)/6)));
    Bxs     = 1/(1+exp((svolt-35)/15));
    Bxs_beta     = 1/(1+exp((svolt-35+IKs_beta_shift)/15));
    %Tau_Xs  = (Axs*Bxs+80);
    beta_Tau_Xs = ((-1.75)*Axs_beta*Bxs_beta+80)*ISO;
    Tau_Xs = beta_Tau_Xs+2.75*Axs*Bxs+80; % beta adrenergic version

        
    R_Inf = 1/(1+exp((20-svolt)/6));
    Tau_R = 9.5*exp(-((svolt+40)^2)/1800)+0.8;
    switch CellType
        case 'Epicardial'
            S_Inf = 1/(1+exp((svolt+20)/5));
            Tau_S = 85*exp(-((svolt+45)^2)/320)+5/(1+exp((svolt-20)/5))+3;
        case 'Mid-myocardial'
            S_Inf = 1/(1+exp((svolt+20)/5));
            Tau_S = 85*exp(-((svolt+45)^2)/320)+5/(1+exp((svolt-20)/5))+3;
        case 'Endocardial'
            S_Inf = 1/(1+exp((svolt+28)/5));
            Tau_S = 1000*exp(-((svolt+67)^2)/1000)+8;
    end
    
    D_Inf = 1/(1+exp((-8-svolt-beta_d_shift_ICal)/7.5));
    Ad    = 1.4/(1+exp((-35-svolt)/13))+0.25;
    Bd    = 1.4/(1+exp((svolt+5)/5));
    Cd    = 1/(1+exp((50-svolt)/20));
    Tau_D = Ad*Bd*Cd;
    
    F_Inf = 1/(1+exp((svolt+20+beta_f_shift_ICal)/7));
    Af    = 1102.5*exp(-((svolt+27)^2)/225);
    Bf    = 200/(1+exp((20.6-svolt)/10));%  % 200/(1+exp((13-svolt)/10))
    Cf    = (180/(1+exp((svolt+30)/10)))+20;
    Tau_F = (Af+Bf+Cf);
    
    F2_Inf = 0.67/(1+exp((svolt+35+beta_f_shift_ICal)/7))+0.33;
    Af2    = 600*exp(-((svolt+25)^2)/170);
    Bf2    = 0; %31/(1+exp((25-svolt)/10));
    Cf2    = 16/(1+exp((svolt+30)/10));
    Tau_F2 = (Af2+Bf2+Cf2);
    
    FCaSS_Inf = 0.6/(1+400*CaSS*CaSS)+0.4;
    Tau_FCaSS = 80/(1+400*CaSS*CaSS)+2;
    
    % Update gates
    sm     = M_Inf-(M_Inf-sm)*exp(-HT/Tau_M);
    %sh     = H_Inf-(H_Inf-sh)*exp(-HT/Tau_H);
    sj     = J_Inf-(J_Inf-sj)*exp(-HT/Tau_J);
    shs    = H_Inf-(H_Inf-shs)*exp(-HT/Tau_H_s);
    shf    = H_Inf-(H_Inf-shf)*exp(-HT/Tau_H_f);
    smL    = ML_Inf-(ML_Inf-smL)*exp(-HT/Tau_ML);
    shL    = HL_Inf-(HL_Inf-shL)*exp(-HT/Tau_HL);
    %sxr1   = Xr1_Inf-(Xr1_Inf-sxr1)*exp(-HT/Tau_Xr1);
    %sxr2   = Xr2_Inf-(Xr2_Inf-sxr2)*exp(-HT/Tau_Xr2);
    xrs     = xrss-(xrss-xrs)*exp(-HT/txrs);
    xrf     = xrss-(xrss-xrs)*exp(-HT/txrf);
    xk1     = xk1ss-(xk1ss-xk1)*exp(-HT/txk1);
    sxs    = (Xs_Inf-(Xs_Inf-sxs)*exp(-HT/Tau_Xs));
    ss     = S_Inf-(S_Inf-ss)*exp(-HT/Tau_S);
    sr     = R_Inf-(R_Inf-sr)*exp(-HT/Tau_R);
    sd     = D_Inf-(D_Inf-sd)*exp(-HT/Tau_D);
    sf     = (F_Inf-(F_Inf-sf)*exp(-HT/Tau_F));%*0.99;
    sf2    = (F2_Inf-(F2_Inf-sf2)*exp(-HT/Tau_F2));%*0.99;
    sfcass = FCaSS_Inf-(FCaSS_Inf-sfcass)*exp(-HT/Tau_FCaSS);%*0.99;
    
    % Update voltage
    svolt = svolt - HT*sItot;
    
        Vm(savecount) = svolt;
        SINaCa(savecount) = INaCa;
        SIKr(savecount) = IKr;
        SIKs(savecount) = IKs;
        SIK1(savecount) = IK1;
        SIto(savecount) = Ito;
        SINa(savecount) = INa;
        SIbNa(savecount) = IbNa;
        SICaL(savecount) = ICaL;
        SIbCa(savecount) = IbCa;
        SINaK(savecount) = INaK;
        SIpCa(savecount) = IpCa;
        SIpK(savecount) = IpK;
        SKi(savecount) = Ki;
        SCai(savecount) = Cai;
        SNai(savecount) = Nai;
        SINaL(savecount) = INaL;
        SItot(savecount) = sItot;
        TimeArray(savecount) = time;
        savecount = savecount + 1;

end % while time

%--- THE MAIN LOOP ENDS HERE ---

savecount = savecount-1;
t = TimeArray(1:savecount);
V = Vm(1:savecount);

%% plot (SHOULD BE DONE IN THE CALLING SCRIPT INSTEAD)
% figure
% 
% plot(V(savecount-55000:savecount-20000), SIK1(savecount-55000:savecount-20000))
% 
% figure;
% subplot(2,2,1); plot(t,V,'k'); title('Action potential');
% subplot(2,2,2); hold on; title('Na+ currents'); % Sodium currents
% plot(t,SINa(1:savecount),'k'); plot(t,SIbNa(1:savecount),'b');
% legend('INa','INaL');
% subplot(2,2,3); hold on; title('K+ currents'); % Potassium currents
% plot(t,SIKr(1:savecount),'k'); plot(t,SIKs(1:savecount),'b'); plot(t,SIK1(1:savecount),'r');
% plot(t,SIto(1:savecount),'g'); plot(t,SIpK(1:savecount),'m');
% legend('IKr','IKs','IK1','Ito','IpK');
% subplot(2,2,4); hold on; title('Ca2+ currents and NaK pump'); % Calcium etc
% plot(t,SICaL(1:savecount),'k'); plot(t,SINaCa(1:savecount),'b'); plot(t,SIpCa(1:savecount),'r');
% plot(t,SIbCa(1:savecount),'g'); plot(t,SINaK(1:savecount),'m');
% legend('ICa','INaCa','IpCa','IbCa','INaK');
% figure;
% plot(t,SItot(1:savecount),'k'); hold on;
% plot(t(1:savecount-1),diff(-Vm(1:savecount)),'b'); hold off;

% figure;
% plot(t,SKi(1:savecount),'k'); hold on;
% plot(t,SCai(1:savecount),'r');
% plot(t,SNai(1:savecount),'b'); hold off;


end % function TenTusscher2

%% EOF