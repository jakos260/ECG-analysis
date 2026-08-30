function result = get_init_rep_combined(GEOM, init_dep, direction_vector, velocity)
    % GET_INIT_REP_COMBINED - Initializes repolarization times with a linear
    % gradient and assigns specific TMP LUTs (EPI/ENDO). Adds a +50 ms offset 
    % to repolarization time for ENDO nodes.

    % 1. Input validation & default values
    if nargin < 4 || isempty(velocity)
        velocity = 1.0; 
    end
    if nargin < 3 || isempty(direction_vector) || numel(direction_vector) ~= 3
        error('A valid 1x3 or 3x1 direction_vector must be provided.');
    end
    if nargin < 2 || isempty(init_dep)
        warning('Depolarization times (init_dep) not provided. result.init_dep will not be set.');
        init_dep = [];
    elseif length(init_dep) ~= length(GEOM.VER)
        error('Length of init_dep must match the number of nodes in GEOM.VER.');
    end
    
    % 2. Calculate the spatial projection for linear gradient
    u = direction_vector(:);
    u = u / norm(u);
    
    projections = GEOM.VER * u;
    distances = projections - min(projections);
    
    % Calculate intermediate dep for rep calculations (1.0 m/s = 1.0 mm/ms)
    dep = distances / velocity; 
    
    % 3. Calculate repolarization times using BSM RMS peaks
    rms_bsm = rms(GEOM.BSM, 1);
    [~, r_time] = max(rms_bsm);
    
    search_start = min(length(rms_bsm), r_time + 100);
    [~, t_time_offset] = max(rms_bsm(search_start:end));
    t_time = search_start - 1 + t_time_offset;
    
    rt = t_time - r_time;
    alpha = 0.8;
    
    % Base init_rep with linear gradient
    init_rep = dep + (rt * alpha);
    
    % 4. Assign potentials (LUTs) and apply ENDO offset
    num_nodes = length(GEOM.VER);
    assigned_LUT = cell(num_nodes, 1);
    
    for i = 1:num_nodes
        if GEOM.endoVER(i) == 0
            assigned_LUT{i} = GEOM.LUT_EPI;
        elseif GEOM.endoVER(i) == 1 || GEOM.endoVER(i) == 2
            assigned_LUT{i} = GEOM.LUT_ENDO;
            % Apply the repolarization offset for ENDO nodes
            init_rep(i) = init_rep(i) + 20;
        else
            assigned_LUT{i} = GEOM.LUT_EPI; % Fallback
        end
    end
    
    % 5. Return structured results
    result = struct();
    if ~isempty(init_dep)
        result.init_dep = init_dep(:);
    end
    result.init_rep = init_rep(:);
    result.assigned_LUT = assigned_LUT;
    
    % Minimal logging
    fprintf('COMBINED SCAN | Vel: %.2f m/s | REP min/max: %4.1f / %4.1f ms\n', ...
            velocity, min(init_rep), max(init_rep));
end