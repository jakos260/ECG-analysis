path = "C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\";
N = size(dir(path), 1) - 2; % remove '.' and '..'

subject_num = 15;
if true 
% for subject_num = 1:N
    subject_name = sprintf('subject_%03d', subject_num);
    % printProgress(subject_num, N, subject_name);

    SubjectPath = fullfile(path, subject_name);
    if isfolder(SubjectPath)
        DATA = readSubject(SubjectPath);
        isSubject = true;
    end

    RESULT = process_subject(DATA);
end


function RESULT = process_subject(DATA)
    GEOM.type = 'ventricles';
    lead_system = 'xstandard12lead'; %'x65Nijmegen' 'x12plus3leads' 'xstandard12lead';
    ref_signal = sprintf('signal_%d', 4); % 3 for LV pacing or 4 for RV pacing
    
    GEOM.BSM        = DATA.SIGNALS.ECG_BEATS.VMEDIANECG.(ref_signal);

    GEOM.VER        = DATA.VENTR.geom.VER;
    GEOM.ITRI       = DATA.VENTR.geom.ITRI;
    GEOM.AMA        = DATA.VENTR.AMA.(lead_system);
    GEOM.DIST       = DATA.VENTR.DIST3D;
    GEOM.DISTsurf   = DATA.VENTR.DISTsurf;
    GEOM.DIST2W     = DATA.VENTR.DISTanis;
    GEOM.ADJ2W      = DATA.VENTR.ADJanis;   % ventricles.adjanis; 
    GEOM.neigh      = DATA.VENTR.ADJneigh;  % ventricle.adjneigh
    GEOM.ADJ        = DATA.VENTR.ADJ3D;     % ventricle.adj3d
    GEOM.LAY        = loadmat('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseArno\BEM\inverse\mla\prague99.mla');
    GEOM.RegionIdx  = DATA.GEOM.ventr.segments;
    
    num_nodes           = length(GEOM.VER);
    GEOM.TVER           = zeros(size(GEOM.BSM, 1), 3);
    GEOM.purkinjever    = zeros(num_nodes, 1);
    GEOM.Lpurkinjever   = zeros(num_nodes, 1);
    GEOM.Rpurkinjever   = zeros(num_nodes, 1);
    GEOM.Rfreewallver   = zeros(num_nodes, 1);
    GEOM.endoVER        = zeros(num_nodes, 1);
    
    
    GEOM.pS = [];
    GEOM.anisotropyRatio = 0.5; % value only to save in logs
    
    % GEOM.specs(2) - początek analizowanego sygnału (np. początek QRS)
    % GEOM.specs(3) - koniec QRS (na tej podstawie liczone jest qrsduration)
    % GEOM.specs(4) - używane do estymacji repolaryzacji
    % GEOM.specs(5) - koniec analizowanego sygnału (koniec załamka T)
    GEOM.specs = [0, 1, 100, 200, 500];
    GEOM.SPECS = struct();
    GEOM.SPECS.onsetqrs = GEOM.specs(2);
    GEOM.SPECS.onsetp = 0; % used only in atria mode
    GEOM.SPECS.qrstduration = GEOM.specs(5) - GEOM.specs(2) + 1;
    GEOM.SPECS.time_Jpoint = GEOM.specs(3);
    GEOM.SPECS.time_Vstim = 1; % used only in stim mode
    GEOM.SPECS.time_apexT = 275;
    GEOM.SPECS.depSlope = 1.0;
    GEOM.SPECS.repCorrection = 0;
    GEOM.SPECS.endtwave = 400;
    GEOM.SPECS.useCumsum = false;
    GEOM.SPECS.initialSlope = 1;
    GEOM.SPECS.plateauslope = 0.02;
    GEOM.SPECS.repslope = 0.045;
    GEOM.SPECS.scaleAmpl = 15;
    load('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseMy\getApLut\ApLut_niceApd_CT1.mat');
    GEOM.LUT = LUT;
    
    % __________initial dep and rep__________
    velocity = 1.0; % m/s
    init = simplyfocusscan(GEOM, velocity);
    init_dep = init.dep;
    
    % init = lineargradientscan(GEOM, [0, 1, -0.2], velocity); % dep vector
    init = lineargradientscan(GEOM, [-1, 0, 0], 0.5);
    init_rep = init.rep;
    
    % _________inv_model setup____________
    mudep = 1e-3;
    murep = 2e-4;
    
    param_configs(1).name = 'Depolarization';
    param_configs(1).initial_value = init_dep;
    param_configs(1).mu = mudep;
    param_configs(1).bounds = [0, 150];
    
    param_configs(2).name = 'Repolarization';
    param_configs(2).initial_value = init_rep;
    param_configs(2).mu = murep;
    param_configs(2).bounds = [150, 400];
    
    settings = struct();
    settings.maxiter = 4;
    settings.plot_en = true;
    settings.lpass = 1;
    settings.mode = 4;
    settings.useWeighedRd = 0;
    settings.MINRD = 0.18;
    settings.stopcrit = 1e-4;
    
    opt_func = @gettres_v_nparams;
    
    inv_model   = INVERSE_MODEL(GEOM, param_configs, settings, opt_func);
    inv_results = inv_model.run_optimization();

    RESULT = struct();
    RESULT.inv_resaults = inv_results;
end