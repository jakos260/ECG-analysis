function [t,V]=TenTusscher2(Stim_I, Stim_T, Stim_Int, CT)
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
%     Stim_I   : Array with the stimulus strengths of monophasic square
%                stimulus pulses (A/F). One entry per pulse.
%     Stim_T   : Array with the respective pulse durations (ms)
%     Stim_Int : Array with interpulse intervals (ms). The interpulse
%                interval is the interval after the start of a pulse to the
%                start of the next pulse or to the end of the simulation.
%     CT       : Cell type: 1=epicardial, 2=mid-myocardial, 3=endocardial
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

Vmaxup = 0.006375;% Maximal Iup conductance (mM/ms)
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
GKr = 0.153;      % Maximal IKr conductance (nS/pF)
switch CellType
    case 'Epicardial',     GKs = 0.392; % Maximal IKs conductance (nS/pF)
    case 'Mid-myocardial', GKs = 0.098; % "
    case 'Endocardial',    GKs = 0.392; % "
end
GK1 = 5.405;        % Maximal IK1 conductance (nS/pF) 
switch CellType
    case 'Epicardial',     Gto = 0.294; % Maximal Ito conductance (nS/pF)
    case 'Mid-myocardial', Gto = 0.294; % "
    case 'Endocardial',    Gto = 0.073; % "
end
GNa   = 14.838;     % Maximal INa conductance (nS/pF)
GbNa  = 0.00029;    % Maximal background-INa conductance (nS/pF)
KmK   = 1.0;
KmNa  = 40.0;
knak  = 2.724;
pKNa  = 0.03;
GCaL  = 0.0000398;
GbCa  = 0.000592;
knaca = 1000;
KmNai = 87.5;
KmCa  = 1.38;
ksat  = 0.1;
n     = 0.35;       % Voltage dependence parameter of INaCa
GpCa  = 0.1238;
KpCa  = 0.0005;
GpK   = 0.0146;

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
sm     = 0.0;
sh     = 0.75;
sj     = 0.75;
sxr1   = 0.0;
sxr2   = 1.0;
sxs    = 0.0;
sr     = 0.0;
ss     = 1.0;
sd     = 0.0;
sf     = 1.0;
sf2    = 1.0;
sfcass = 1.0;
sRR    = 1.0;
sOO    = 0.0;

% End of Electrophysiological parameters

%% Integration parameters

HT = 0.02;                       % timestep (ms)
STOPTIME = 1000 + sum(Stim_Int); % Duration of the simulation (ms)

% End of Integration parameters

LEN = ceil(STOPTIME/HT)+1;
Vm = zeros(LEN,1);         % Membrane potential (mV)
SINaCa = zeros(LEN,1);     % INaCa to be saved
SIKr = zeros(LEN,1);
SIKs = zeros(LEN,1);
SIK1 = zeros(LEN,1);
SIto = zeros(LEN,1);
SINa = zeros(LEN,1);
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
    AK1 = 0.1/(1+exp(0.06*(svolt-EK-200)));
    BK1 = (3*exp(0.0002*(svolt-EK+100)) + exp(0.1*(svolt-EK-10))) / ...
           (1+exp(-0.5*(svolt-EK)));
    rec_iK1  = AK1/(AK1+BK1);
    rec_iNaK = 1 / ...
                (1+0.1245*exp(-0.1*svolt/RTonF)+0.0353*exp(-svolt/RTonF));
    rec_ipK = 1 / (1+exp((25-svolt)/5.98));
    
    % Calculate ion currents
    INa  = GNa*(sm^3)*sh*sj*(svolt-ENa);
    svolt_etc = 2*(svolt-15)/RTonF;
    ICaL = GCaL*sd*sf*sf2*sfcass*2*F*svolt_etc * ...
           (0.25*exp(svolt_etc)*CaSS-Cao) / (exp(svolt_etc)-1);
    Ito = Gto*sr*ss*(svolt-EK);
    IKr = GKr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-EK);
    IKs = GKs*sxs*sxs*(svolt-EKs);
    IK1 = GK1*rec_iK1*(svolt-EK);
    INaCa = knaca*(1/((KmNai^3)+(Nao^3)))*(1/(KmCa+Cao))* ...
            (1/(1+ksat*exp((n-1)* svolt/RTonF)))* ...
            (exp(n*svolt/RTonF)*(Nai^3)*Cao-exp((n-1)*svolt/RTonF)* ...
            (Nao^3)*Cai*2.5);
    INaK = knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
%    INaK = 0.5*knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa = GpCa*Cai/(KpCa+Cai);
    IpK  = GpK*rec_ipK*(svolt-EK);
    IbNa = GbNa*(svolt-ENa);
    IbCa = GbCa*(svolt-ECa);
    
    % Total current is sum of all the above
    sItot = IKr+IKs+IK1+Ito+INa+IbNa+ICaL+IbCa+INaK+INaCa+IpCa+IpK+Istim(savecount);
    
    % Update concentrations
    kCaSR = maxsr-(maxsr-minsr)/(1+(EC/CaSR)^2);
    k1    = k1_/kCaSR;
    k2    = k2_*kCaSR;
    dRR   = k4*(1-sRR)-k2*CaSS*sRR;
    sRR   = sRR+HT*dRR;
    sOO   = k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);
    
    Irel  = Vrel*sOO*(CaSR-CaSS);
    Ileak = Vleak*(CaSR-Cai);
    Iup   = Vmaxup/(1+(Kup/Cai)^2);
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
    
    CaBuf = Bufc*Cai/(Kbufc+Cai);
    dCai  = HT*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*Capacitance)- ...
             (Iup-Ileak)*(Vsr/Vc)+Ixfer);
    bc    = Bufc-CaBuf-dCai-Cai+Kbufc;
    cc    = Kbufc*(CaBuf+dCai+Cai);
    Cai   = (sqrt(bc*bc+4*cc)-bc)/2;
    
    dNai = -(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*Capacitance;
    Nai  = Nai + HT*dNai;
    
    dKi = -(Istim(savecount)+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*Capacitance;
    Ki  = Ki + HT*dKi;
        
    % Calculate steady state values and time constants (gating parameters)
    AM    = 1/(1+exp((-60-svolt)/5));
    BM    = 0.1/(1+exp((svolt+35)/5))+0.1/(1+exp((svolt-50)/200));
    Tau_M = AM*BM;
    M_Inf = 1/((1+exp((-56.86-svolt)/9.03))^2);
    if svolt>-40
        AH_1  = 0;
        BH_1  = 0.77/(0.13*(1+exp(-(svolt+10.66)/11.1)));
        Tau_H =1/(AH_1+BH_1);
    else
        AH_2  = (0.057*exp(-(svolt+80)/6.8));
        BH_2  = 2.7*exp(0.079*svolt)+3.1e5*exp(0.3485*svolt);
        Tau_H = 1/(AH_2+BH_2);
    end
    H_Inf = 1/((1+exp(svolt+71.55)/7.43))^2;
    if svolt>-40
        AJ_1  = 0;
        BJ_1  = (0.6*exp((0.057)*svolt)/(1+exp(-0.1*(svolt+32))));
        Tau_J = 1/(AJ_1+BJ_1);
    else
        AJ_2  = (((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)* ...
                 exp(-0.04391*svolt))*(svolt+37.78)/ ...
                 (1+exp(0.311*(svolt+79.23))));
        BJ_2  = (0.02424*exp(-0.01052*svolt)/ ...
                 (1+exp(-0.1378*(svolt+40.14))));
        Tau_J = 1/(AJ_2+BJ_2);
    end
    J_Inf = H_Inf;
    
    Xr1_Inf = 1/(1+exp((-26-svolt)/7));
    axr1    = 450/(1+exp((-45-svolt)/10));
    bxr1    = 6/(1+exp((svolt+30)/11.5));
    Tau_Xr1 = axr1*bxr1;
    
    Xr2_Inf = 1/(1+exp((svolt+88)/24));
    axr2    = 3/(1+exp((-60-svolt)/20));
    bxr2    = 1.12/(1+exp((svolt-60)/20));
    Tau_Xr2 = axr2*bxr2;
    
    Xs_Inf  = 1/(1+exp((-5-svolt)/14));
    Axs     = 1400/(sqrt(1+exp((5-svolt)/6)));
    Bxs     = 1/(1+exp((svolt-35)/15));
    Tau_Xs  = Axs*Bxs+80;
        
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
    
    D_Inf = 1/(1+exp((-8-svolt)/7.5));
    Ad    = 1.4/(1+exp((-35-svolt)/13))+0.25;
    Bd    = 1.4/(1+exp((svolt+5)/5));
    Cd    = 1/(1+exp((50-svolt)/20));
    Tau_D = Ad*Bd*Cd;
    
    F_Inf = 1/(1+exp((svolt+20)/7));
    Af    = 1102.5*exp(-((svolt+27)^2)/225);
    Bf    = 200/(1+exp((13-svolt)/10));
    Cf    = (180/(1+exp((svolt+30)/10)))+20;
    Tau_F = Af+Bf+Cf;
    
    F2_Inf = 0.67/(1+exp((svolt+35)/7))+0.33;
    Af2    = 600*exp(-((svolt+25)^2)/170);
    Bf2    = 31/(1+exp((25-svolt)/10));
    Cf2    = 16/(1+exp((svolt+30)/10));
    Tau_F2 = Af2+Bf2+Cf2;
    
    FCaSS_Inf = 0.6/(1+400*CaSS*CaSS)+0.4;
    Tau_FCaSS = 80/(1+400*CaSS*CaSS)+2;
    
    % Update gates
    sm     = M_Inf-(M_Inf-sm)*exp(-HT/Tau_M);
    sh     = H_Inf-(H_Inf-sh)*exp(-HT/Tau_H);
    sj     = J_Inf-(J_Inf-sj)*exp(-HT/Tau_J);
    sxr1   = Xr1_Inf-(Xr1_Inf-sxr1)*exp(-HT/Tau_Xr1);
    sxr2   = Xr2_Inf-(Xr2_Inf-sxr2)*exp(-HT/Tau_Xr2);
    sxs    = Xs_Inf-(Xs_Inf-sxs)*exp(-HT/Tau_Xs);
    ss     = S_Inf-(S_Inf-ss)*exp(-HT/Tau_S);
    sr     = R_Inf-(R_Inf-sr)*exp(-HT/Tau_R);
    sd     = D_Inf-(D_Inf-sd)*exp(-HT/Tau_D);
    sf     = F_Inf-(F_Inf-sf)*exp(-HT/Tau_F);
    sf2    = F2_Inf-(F2_Inf-sf2)*exp(-HT/Tau_F2);
    sfcass = FCaSS_Inf-(FCaSS_Inf-sfcass)*exp(-HT/Tau_FCaSS);
    
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
        SItot(savecount) = sItot;
        TimeArray(savecount) = time;
        savecount = savecount + 1;

end % while time

%--- THE MAIN LOOP ENDS HERE ---

savecount = savecount-1;
t = TimeArray(1:savecount);
V = Vm(1:savecount);

%% plot (SHOULD BE DONE IN THE CALLING SCRIPT INSTEAD)

figure;
subplot(2,2,1); plot(t,V,'k'); title('Action potential');
subplot(2,2,2); hold on; title('Na+ currents'); % Sodium currents
plot(t,SINa(1:savecount),'k'); plot(t,SIbNa(1:savecount),'b');
legend('INa','INaL');
subplot(2,2,3); hold on; title('K+ currents'); % Potassium currents
plot(t,SIKr(1:savecount),'k'); plot(t,SIKs(1:savecount),'b'); plot(t,SIK1(1:savecount),'r');
plot(t,SIto(1:savecount),'g'); plot(t,SIpK(1:savecount),'m');
legend('IKr','IKs','IK1','Ito','IpK');
subplot(2,2,4); hold on; title('Ca2+ currents and NaK pump'); % Calcium etc
plot(t,SICaL(1:savecount),'k'); plot(t,SINaCa(1:savecount),'b'); plot(t,SIpCa(1:savecount),'r');
plot(t,SIbCa(1:savecount),'g'); plot(t,SINaK(1:savecount),'m');
legend('ICa','INaCa','IpCa','IbCa','INaK');
% figure;
% plot(t,SItot(1:savecount),'k'); hold on;
% plot(t(1:savecount-1),diff(-Vm(1:savecount)),'b'); hold off;

% figure;
% plot(t,SKi(1:savecount),'k'); hold on;
% plot(t,SCai(1:savecount),'r');
% plot(t,SNai(1:savecount),'b'); hold off;


end % function TenTusscher2

%% EOF