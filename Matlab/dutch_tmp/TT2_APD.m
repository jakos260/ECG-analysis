% TT2_APD  To test the Ten Tusscher Mark 2 model for APD dependency of the
%          stimulus interval
% For rate below 60 bpm the APD is almost constant, which is a flaw in
% the model. For rate > 60 bpm there is a clear change in APD with rate.
%%
I0  = -24;  % Stimulus strength (A/F)
T0  = 1;    % Pulse duration (ms)
RR0 = 1000; % Standard RR interval (ms)

%% Create pulse train

% Normalized stimulus current for each pulse:
I = [1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% Normalized pulse duration for each pulse:
T = [1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% Normalized interpulse interval after start of each pulse
RR = [0.1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% ( RR(1) is the starttime of the first pulse)

%% Run Ten Tusscher 2 model with the pulse train
disp('TT2_APD running');
for CellType = 1:3 % epicardial, M-cell and endocardial
    [t, VAP] = TenTusscher2(I0*I, T0*T, RR0*RR, CellType);
    n1 = length(t)-55000;
    n2 = length(t)-25000;
    t_singleAP = t(n1:n2)-t(n1);
    V_singleAP = VAP(n1:n2);
    switch CellType
        case 1 % Epicardial
            figure;
            plot(t_singleAP, V_singleAP,'k');
            title('Epicardial');
            xlabel('t (ms)');
            ylabel('Vmembrane (mV)');
            save("AP_epicardial.txt","t_singleAP","V_singleAP","-ascii");
        case 2 % Mid-myocardial
            figure;
            plot(t_singleAP, V_singleAP,'k');
            title('Mid-myocardial');
            xlabel('t (ms)');
            ylabel('Vmembrane (mV)');
            save("AP_midmyocardial.txt","t_singleAP","V_singleAP","-ascii");
        case 3 % Endocardial
            figure;
            plot(t_singleAP, V_singleAP,'k');
            title('Endocardial');
            xlabel('t (ms)');
            ylabel('Vmembrane (mV)');
            save("AP_endocardial.txt","t_singleAP","V_singleAP","-ascii");
    end
end

%% EOF