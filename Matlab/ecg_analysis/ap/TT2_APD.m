% TT2_APD  To test the Ten Tusscher Mark 2 model for APD dependency of the
%          stimulus interval
% For rate below 60 bpm the APD is almost constant, which is a flaw in
% the model. For rate > 60 bpm there is a clear change in APD with rate.
%
APs = zeros();
for pcl = 1
    I0  = -30;  % Stimulus strength (A/F)
    T0  = 1;    % Pulse duration (ms)
    ISO = [0.1];  % Ratio of beta-adrenergic activation
    RR0 = 1000*(60/(120*ISO(pcl)+60)); % Standard RR interval (ms) 1000 = 60 bpm
    % ISO = 0 -> RR0 = 1000 = 60 bpm. ISO = 1 -> RR0 = 333.3 = 180 bpm
    
    I = ones(1,51);
    % Normalized stimulus current for each pulse:
    %I = [1 1 1 1 1 1 1 1 1 1 1 1 1 ];
    % Normalized pulse duration for each pulse:
    T = ones(1,51);
    %T = [1 1 1 1 1 1 1 1 1 1 1 1 1 ];
    % Normalized interpulse interval after start of each pulse
    RR = [0.1 ones(1,50)];
    %RR = [0.1 1 1 1 1 1 1 1 1 1 1 1 1 ];
    % ( RR(1) is the starttime of the first pulse)
    
    % Run Ten Tusscher 2 model with the pulse train
    disp('TT2_APD running');
    for CellType = 1:3 % epicardial, M-cell and endocardial
        [t, VAP] = TenTusscher2(I0*I, T0*T, ISO(pcl), RR0*RR, CellType);
        n1 = length(t)-50*1000; %50150;
        n2 = length(t)-50*500; %23000;
        t_singleAP = t(n1:n2)-t(n1);
        V_singleAP = VAP(n1:n2);
        index_of_stim = cumsum(RR0*RR);
        last_stim_time = index_of_stim(end);
        t_for_last_stim = last_stim_time*50-n1;
        switch CellType
            case 1 % Epicardial
                figure;
                plot(t_singleAP, V_singleAP,'k');
                title('Epicardial');
                xlabel('t (ms)');
                ylabel('Vmembrane (mV)');
                % Calculate APD90
                min_single_AP = min(V_singleAP);
                Positive_AP = (-min_single_AP)+V_singleAP;
                idx = find(Positive_AP(200:end) <= max(Positive_AP(200:end)*0.1));
                APD90_index = idx(1)+200;
                APD90_index_real = (APD90_index-t_for_last_stim)/50;
                fprintf('APD90 for Epicardial cell = %.4f ms\n', APD90_index_real);
                % corrected APD
                cAPD = (APD90_index_real/1000)/((RR0/1000)^(1/3));
                fprintf('corrected APD90 for Epicardial cell = %.4f ms\n', cAPD*1000);
            case 2 % Mid-myocardial
                figure;
                plot(t_singleAP, V_singleAP,'k');
                title('Mid-myocardial');
                % xlim([0 1000])  % window limits
                % ylim([-100 60]) % window limits
                xlabel('t (ms)');
                ylabel('Vmembrane (mV)');
                % Calculate APD90
                min_single_AP = min(V_singleAP);
                Positive_AP = (-min_single_AP)+V_singleAP;
                idx = find(Positive_AP(200:end) <= max(Positive_AP(200:end)*0.1));
                APD90_index = idx(1)+200;
                APD90_index_real = (APD90_index-t_for_last_stim)/50;
                fprintf('APD90 for Midmyocardial cell = %.4f ms\n', APD90_index_real);
                % corrected APD
                cAPD = (APD90_index_real/1000)/((RR0/1000)^(1/3));
                fprintf('corrected APD90 for Epicardial cell = %.4f ms\n', cAPD*1000);
            case 3 % Endocardial
                figure;
                plot(t_singleAP, V_singleAP,'k');
                title('Endocardial');
                xlabel('t (ms)');
                ylabel('Vmembrane (mV)');
                % Calculate APD90
                min_single_AP = min(V_singleAP);
                Positive_AP = (-min_single_AP)+V_singleAP;
                idx = find(Positive_AP(200:end) <= max(Positive_AP(200:end)*0.1));
                APD90_index = idx(1)+200;
                APD90_index_real = (APD90_index-t_for_last_stim)/50;
                fprintf('APD90 for endocardial cell = %.4f ms\n', APD90_index_real);
                % corrected APD
                cAPD = (APD90_index_real/1000)/((RR0/1000)^(1/3));
                fprintf('corrected APD90 for Epicardial cell = %.4f ms\n', cAPD*1000);
        end
    end
    APs(pcl, 1:length(V_singleAP)) = V_singleAP;
end
