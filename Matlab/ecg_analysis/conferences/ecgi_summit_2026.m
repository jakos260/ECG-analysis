path = "C:\Users\Admin\Documents\Projects\ecg_project\Scripts\data\Dataset\";
N = size(dir(path), 1) - 2; % remove '.' and '..'

subject_num = 15;
if true 
% for subject_num = 1:N
    subject_name = sprintf('subject_%03d', subject_num);
    % printProgress(subject_num, N, subject_name);

    SubjectPath = fullfile(path, subject_name);
    DATA = readSubject(SubjectPath);

    RESULT = process_subject(DATA);
end

%%
show_result(DATA, RESULT);

%%
function show_result(DATA, RESULT)
    qtripy_en = true;
    ecg_en = true;

    if ecg_en
        signal_name     = RESULT.signal_name;
        LAY             = loadmat('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseArno\BEM\inverse\mla\standard12.mla');
        onset_qrs       = DATA.SIGNALS.INFO.(signal_name).MORPHOLOGY.ALL.SPECS.onsetqrs;
        end_twave       = DATA.SIGNALS.INFO.(signal_name).MORPHOLOGY.ALL.SPECS.endtwave;

        signal_ref      = DATA.SIGNALS.ECG_BEATS.VMEDIANECG.(signal_name)(:, onset_qrs:end_twave);
        signal_inv      = RESULT.inv_signal;
        signal_map      = RESULT.map_signal;

        figure(1); clf; hold on;
        set(gcf, 'Name', 'Result');
        sigplot(signal_ref, 'Black = Measured, Red = Simulated, Blue = Mapper', LAY, 1.5, 'k', 1);
        sigplot(lowpassma(signal_inv, 5), '', LAY, 1.5, 'r', 1);
        sigplot(lowpassma(signal_map, 5), '', LAY, 1.5, 'b', 1);

        fprintf("Result RD new_inv: %.3f mapper:%.3f\n", RelativeDistance(signal_ref, signal_inv), RelativeDistance(signal_ref, signal_map));
    end

    if qtripy_en
        transparency = 0.0;
        gradient_bins = 25;
        q = initQtripy();
        
        q.reset();
        q.disable_debounce();
        q.set_panels_number(1,3);
        q.background_color("white");
        q.color_range(200, 400);
        q.gradient_bins(gradient_bins);
        q.text(RESULT.signal_name, [0.3, 0.95]);
        
        % mapper
        q.set_active_panel(1, 1);
        q.surface(DATA.GEOM.ventr.VER, DATA.GEOM.ventr.ITRI);
        q.values(RESULT.map_rep);
        q.text(sprintf("mapper_rep"), [0.1, 0.9]);
        
        % inv
        q.set_active_panel(1, 2);
        q.surface(DATA.GEOM.ventr.VER, DATA.GEOM.ventr.ITRI);
        q.values(RESULT.inv_rep);
        q.text(sprintf("inv(ap)_rep"), [0.1, 0.6]);

        % init_rep
        q.set_active_panel(1, 3);
        q.surface(DATA.GEOM.ventr.VER, DATA.GEOM.ventr.ITRI);
        q.values(RESULT.init_rep);
        q.text(sprintf("init_rep"), [0.1, 0.3]);

    end
end


function RESULT = process_subject(DATA)
    GEOM.type = 'ventricles';
    lead_system = 'xstandard12lead'; %'x65Nijmegen' 'x12plus3leads' 'xstandard12lead';
    ref_signal = sprintf('signal_%d', 4); % 3 for LV pacing or 4 for RV pacing
    
    onset_qrs       = DATA.SIGNALS.INFO.(ref_signal).MORPHOLOGY.ALL.SPECS.onsetqrs;
    end_twave       = DATA.SIGNALS.INFO.(ref_signal).MORPHOLOGY.ALL.SPECS.endtwave;
    duration_qrs    = DATA.SIGNALS.INFO.(ref_signal).MORPHOLOGY.ALL.SPECS.qrsduration;
    twave_peak      = DATA.SIGNALS.INFO.(ref_signal).MORPHOLOGY.ALL.SPECS.time_apexT - onset_qrs;
    GEOM.BSM        = baselinecor(DATA.SIGNALS.ECG_BEATS.VMEDIANECG.(ref_signal)(:, onset_qrs:end_twave));

    GEOM.VER        = DATA.VENTR.geom.VER;
    GEOM.ITRI       = DATA.VENTR.geom.ITRI;
    GEOM.AMA        = DATA.VENTR.AMA.(lead_system);
    GEOM.DIST       = DATA.VENTR.DIST3D;
    GEOM.DISTsurf   = DATA.VENTR.DISTsurf;
    GEOM.DIST2W     = DATA.VENTR.DISTanis;
    GEOM.ADJ2W      = DATA.VENTR.ADJanis;   % ventricles.adjanis; 
    GEOM.neigh      = DATA.VENTR.ADJneigh;  % ventricle.adjneigh
    GEOM.ADJ        = DATA.VENTR.ADJ3D;     % ventricle.adj3d
    GEOM.LAY        = loadmat('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseArno\BEM\inverse\mla\standard12.mla'); % prague99.mla
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
    GEOM.SPECS.signalLength = size(GEOM.BSM, 2);
    GEOM.SPECS.useCumsum = false;
    GEOM.SPECS.initialSlope = 1;
    GEOM.SPECS.plateauslope = 0.02;
    GEOM.SPECS.repslope = 0.045;
    GEOM.SPECS.scaleAmpl = 15;

    load('C:\Users\Admin\Documents\Projects\ecg_project\Scripts\Matlab\ecg_analysis\inverseMy\getApLut\ApLut_niceApd.mat');
    lower_lim = twave_peak - duration_qrs * 0.25;
    upper_lim = twave_peak + duration_qrs * 0.25;
    GEOM.LUT = LUT([LUT.APD] > lower_lim & [LUT.APD] < upper_lim);
    
    % __________initial dep and rep__________
    velocity = 1.0; % m/s
    init = simplyfocusscan(GEOM, velocity);
    init_dep = DATA.MAPPER.x12ECG.(ref_signal).BEATS.ALL.RES(1).finaldep;
    
    % init = lineargradientscan(GEOM, [0, 1, -0.2], velocity); % dep vector
    % init = lineargradientscan(GEOM, [-1, 0, 0], 0.5);
    init = lineargradientscan(GEOM, [-1, 0, 0], [lower_lim, upper_lim]);
    init_rep = init.rep;
    
    % _________inv_model setup____________
    murep = 1e-4;
    
    param_configs(1).name = 'Repolarization';
    param_configs(1).initial_value = init_rep;
    param_configs(1).mu = murep;
    param_configs(1).bounds = [150, 400];
    
    settings = struct();
    settings.maxiter = 4;
    settings.plot_en = true;
    settings.lpass = 1;
    settings.mode = 4;
    settings.useWeighedRd = 0;
    settings.MINRD = 0.18;
    settings.stopcrit = 1e-4;
    
    fixed_dep = init_dep;
    opt_func = @(INV, params, notchopt, amplopt) gettres_v_rep(INV, params, fixed_dep, notchopt, amplopt);
    
    inv_model   = INVERSE_MODEL(GEOM, param_configs, settings, opt_func);
    inv_results = inv_model.run_optimization();

    opt_func_mapper = @(INV, params, notchopt, amplopt) gettres_v_mapper(INV, params, fixed_dep, notchopt, amplopt);
    mapper_model    = INVERSE_MODEL(GEOM, param_configs, settings, opt_func_mapper);
    mapper_results  = mapper_model.run_optimization();

    RESULT = struct();
    RESULT.signal_name  = ref_signal;
    RESULT.lead_name    = lead_system;
    RESULT.fixed_dep    = fixed_dep;
    RESULT.map_data_dep = DATA.MAPPER.x12ECG.(ref_signal).BEATS.ALL.RES(1).finaldep;
    RESULT.map_data_rep = DATA.MAPPER.x12ECG.(ref_signal).BEATS.ALL.RES(1).finalrep;
    RESULT.inv_rep      = inv_results.final_params{1};
    RESULT.inv_signal   = inv_results.simulated_ecg;
    RESULT.map_rep      = mapper_results.final_params{1};
    RESULT.map_signal   = mapper_results.simulated_ecg;
    RESULT.init_rep     = init_rep;
end