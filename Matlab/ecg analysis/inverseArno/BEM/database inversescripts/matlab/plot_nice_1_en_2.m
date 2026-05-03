
function plot_nice_1_en_2

%% initialize data and set parameters:
clear all;
close all;

% Define subject and beats:
SUB         = [1 2];
beat_nr     = 1:5;
n_filt      = 0;
[bp]        = '/Users/arnojanssen/Documents/STW/PVCs/'; % main directory NICE analyses

for S = SUB,
    
    subject     = ['NICE00' num2str(S)];
    data_dir    = [bp subject '/results/cluster1/'];
    
    if S == 1,
        foldername   = '20150316_mode1_rd0.25_im0_wrd0_iV0.4_iANIS2.0ANIS2.00mudep1.00e-05/';
        fname_p1     = 'NICE001_NICE001_PO_';
        fname_p2     = '_ventricles_im0_wrd0_iV0.4_iANIS2.0ANIS2.00_20150316.mat';
        beat_array   = {'beat031'; 'beat085'; 'beat106'; 'beat196'; 'beat205';};
        
        lay_ch_n     = [100 9 15 100 23 27 44 45 46];
        wct_n        = [44 45 46];
        
        lay_ch_ecg   = [7 8 9];
        nr_sub_ecg   = [2 2];
        sub_ecg_plot = [1 2 4];
        sub_ecg_name = {'D1'; 'E1'; 'E3/V2';};
        
        % set plot parameters:
        xminmax         = [0 500];
        xminmax_short   = [0 110];
        yminmax         = [-2 2.5];
        yminmax_short   = [-0.5 0.5];
        yminmax_ecg     = [-3 1.5];
        %xminmax_ecg     = [0 200];
        
    elseif S == 2,
        foldername   = '20150316_mode1_rd0.25_im0_wrd0_iV0.4_iANIS2.0ANIS2.00mudep1.00e-05/';
        fname_p1     = 'NICE002_NICE002_BSM_0_59_65ch_20150225T121459_';
        fname_p2     = '_ventricles_im0_wrd0_iV0.4_iANIS2.0ANIS2.00_20150316.mat';
        beat_array   = {'beat002'; 'beat060'; 'beat074'; 'beat077'; 'beat083';};
        lay_ch_n     = [12 18 25 31 40 45 61 62 63];
        wct_n        = [61 62 63];
        
        lay_ch_ecg   = [12 13 14 15 18 19 20 21 25 26];
        nr_sub_ecg   = [4 3];
        sub_ecg_plot = [1 4 7 10 2 5 8 11 6 9];
        sub_ecg_name = {'D3/V1'; 'D4'; 'D5'; 'D6'; 'E3/V2'; 'E4'; 'E5'; 'E6'; 'F2/V3'; 'F3'};
        
        % set plot parameters:
        xminmax         = [0 500];
        xminmax_short   = [0 110];
        yminmax         = [-1.5 1.5];
        yminmax_short   = [-0.5 0.5];
        yminmax_ecg     = [-1.5 1.5];
    end
    
    %% load data per beat and plot 12-lead ECG:
    for beats = beat_nr,
        
        % load data
        load([data_dir foldername fname_p1 beat_array{beats} fname_p2])
        
        SIG     = GEOM.BSM;
        SIG_B   = lowpassma(SIG,5);
        SIG_F   = baselinecor(SIG_B);
        
        % apply notch filter:
        if n_filt == 1, SIG_F = subtracthum_aj(SIG_F); nf_tag = 'NF_'; else nf_tag = ''; end
        
        % WCT:
        wct     = mean(SIG_F(wct_n,:));
        
        %% Constrcut 12 Lead plot:
        LI      = SIG_F(lay_ch_n(7),:) - SIG_F(lay_ch_n(8),:);
        LII     = SIG_F(lay_ch_n(9),:) - SIG_F(lay_ch_n(8),:);
        LIII    = SIG_F(lay_ch_n(9),:) - SIG_F(lay_ch_n(7),:);
        avL     = 1.5*(SIG_F(lay_ch_n(7),:) - wct);
        avR     = 1.5*(SIG_F(lay_ch_n(8),:) - wct);
        avF     = 1.5*(SIG_F(lay_ch_n(9),:) - wct);
        V2      = SIG_F(lay_ch_n(2),:) - wct;
        V3      = SIG_F(lay_ch_n(3),:) - wct;
        V6  	= SIG_F(lay_ch_n(6),:) - wct;
        
        % construct V5 and V6: bad channels are empty:
        if S == 2,
            V1      = SIG_F(lay_ch_n(1),:) - wct;
            V4      = SIG_F(lay_ch_n(4),:) - wct;
            V5    	= NaN(1,length(SIG_F));
            
        elseif S == 1,
            V1    	= NaN(1,length(SIG_F));
            V4    	= NaN(1,length(SIG_F));
            V5    	= SIG_F(lay_ch_n(5),:) - wct;
        end
        
        % construct vector with subplot titles and data:
        plot_lead_title = {'LI'; 'avR'; 'V1'; 'V4'; 'LII'; 'avL'; 'V2'; 'V5'; 'LIII'; 'avF'; 'V3'; 'V6'};
        plot_lead       = [LI; avR; V1; V4; LII; avL; V2; V5; LIII; avF; V3; V6];
        
        % plot figures:
        f1 = figure(1);
        for k = 1:12,
            subplot(3,4,k)
            plot(plot_lead(k,:));
            hold on;
            plot(zeros(1,size(plot_lead(k,:),2)), 'r');
            grid minor
            ylim(yminmax)
            xlim(xminmax)
            title(plot_lead_title{k})
            xlabel('time (ms)')
            ylabel('mV')
        end
        
        f2 = figure(2);
        for k = 1:12,
            subplot(3,4,k)
            plot(plot_lead(k,:));
            hold on;
            plot(zeros(1,size(plot_lead(k,:),2)), 'r');
            grid minor
            ylim(yminmax_short)
            xlim(xminmax_short)
            title(plot_lead_title{k})
            xlabel('time (ms)')
            ylabel('mV')
        end
        
        % save figures:
        saveas(f1, [data_dir foldername 'figures/ Twelve_LEAD_' nf_tag subject '_' beat_array{beats} '.jpg'])
        saveas(f2, [data_dir foldername 'figures/ Twelve_LEAD_' nf_tag subject '_' beat_array{beats} '_short.jpg'])
        
        close all
        
        %% Construct ECG plot with channels near origin:
        
        % construct vector with subplot data:
        SIG_F_zm    = zeromean(SIG_F);
        plot_lead   = SIG_F_zm(lay_ch_ecg,:);
        
        % plot figures:
        f3 = figure(3);
        for k = 1:size(plot_lead,1),
            subplot(nr_sub_ecg(1),nr_sub_ecg(2),sub_ecg_plot(k))
            plot(plot_lead(k,:));
            hold on;
            plot(zeros(1,size(plot_lead(k,:),2)), 'r');
            grid minor
            ylim(yminmax_ecg)
            xlim(xminmax)
            title(sub_ecg_name{k})
            xlabel('time (ms)')
            ylabel('mV')
        end
        
        f4 = figure(4);
        for k = 1:size(plot_lead,1),
            subplot(nr_sub_ecg(1),nr_sub_ecg(2),sub_ecg_plot(k))
            plot(plot_lead(k,:));
            hold on;
            plot(zeros(1,size(plot_lead(k,:),2)), 'r');
            grid minor
            ylim(yminmax_short)
            xlim(xminmax_short)
            title(sub_ecg_name{k})
            xlabel('time (ms)')
            ylabel('mV')
        end
        
        % save figures:
        saveas(f3, [data_dir foldername 'figures/ ECG_ORIGIN_' nf_tag subject '_' beat_array{beats} '.jpg'])
        saveas(f4, [data_dir foldername 'figures/ ECG_ORIGIN_' nf_tag subject '_' beat_array{beats} '_short.jpg'])
        
         close all
        
    end
    
end

function sigout = subtracthum_aj(sigin)

fs = 1000;             %#sampling rate
f0 = 50;                %#notch frequency
fn = fs/2;              %#Nyquist frequency
freqRatio = f0/fn;      %#ratio of notch freq. to Nyquist freq.

notchWidth = 0.1;       %#width of the notch

% Compute zeros
zeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

% Compute poles
poles = (1-notchWidth) * zeros;

b = poly( zeros ); %# Get moving average filter coefficients
a = poly( poles ); %# Get autoregressive filter coefficients

% filter signal x
for k = 1: size(sigin,1),
    sigout(k,:) = filter(b,a,sigin(k,:));
end