% %% initialize data and set parameters:
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
        load('/Users/arnojanssen/Documents/STW/PVCs/NICE001/Biosemi/export/NICE001_PO.mat');
        
    elseif S == 2,
        foldername   = '20150316_mode1_rd0.25_im0_wrd0_iV0.4_iANIS2.0ANIS2.00mudep1.00e-05/';
        fname_p1     = 'NICE002_NICE002_BSM_0_59_65ch_20150225T121459_';
        fname_p2     = '_ventricles_im0_wrd0_iV0.4_iANIS2.0ANIS2.00_20150316.mat';
        beat_array   = {'beat002'; 'beat060'; 'beat074'; 'beat077'; 'beat083';};
        load('/Users/arnojanssen/Documents/STW/PVCs/NICE002/Biosemi/export/NICE002_BSM_0_59_65ch_20150225T121459.mat');
        
    end
    
    for beats = beat_nr,
        
        % load data
        load([data_dir foldername fname_p1 beat_array{beats} fname_p2])
        
        % create heart with deps:
        qtriplot(GEOM.VER,GEOM.ITRI)
        qtriplot(meas.depfinal)
        qtriplot('funcolor tim')
        qtriplot('funscale auto')
        qtriplot('bgdcolor white');
        
        % create thorax with elecs:
        qtriplot(GEOM.TVER,GEOM.TITRI)
        qtriplot('color flesh')
        qtriplot('trans 0.25')
        qtriplot(GEOM.TVER(find(DATA.remove == 0)',:),[],'BSMelecs');
        qtriplot('color black');
        
        %pause
        
        % save figure:
        qtriplot(['png ' data_dir foldername 'figures/Thorax_with_elecs_' subject '_' beat_array{beats} '.png 1000 1000']); pause(0.2);
        
        qtriplot('exit')
        
    end
    
end