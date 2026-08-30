function result = get_init_rep_epiendo(GEOM, init_dep)
    % FIXEDREPOLARIZATIONSCAN - Assigns a rigid repolarization time to all nodes
    % based on the BSM T-wave peak and matches specific TMP LUTs (EPI/ENDO)
    % to the provided depolarization times.

    % Check if depolarization times are provided
    if nargin < 2 || isempty(init_dep)
        error('Depolarization times (init_dep) must be provided.');
    end

    % 1. Calculate repolarization time using BSM RMS peak
    rms_bsm = rms(GEOM.BSM, 1);
    [~, r_time] = max(rms_bsm);
    
    % Find T-wave peak after QRS[cite: 1]
    search_start = min(length(rms_bsm), r_time + 100);
    [~, t_time_offset] = max(rms_bsm(search_start:end));
    t_time = search_start - 1 + t_time_offset;
    
    % 2. Assign the fixed repolarization time to all nodes
    num_nodes = length(GEOM.VER);
    init_rep = repmat(t_time, num_nodes, 1);
    
    % Ensure init_dep has the correct dimensions
    if length(init_dep) ~= num_nodes
        error('Length of init_dep must match the number of nodes in GEOM.VER.');
    end
    
    % 3. Assign potentials (LUTs) based on GEOM.endoVER
    assigned_LUT = cell(num_nodes, 1);
    
    for i = 1:num_nodes
        if GEOM.endoVER(i) == 0
            assigned_LUT{i} = GEOM.LUT_EPI;
        elseif GEOM.endoVER(i) == 1 || GEOM.endoVER(i) == 2
            assigned_LUT{i} = GEOM.LUT_ENDO;
            init_rep(i) = init_rep(i) + 50;
        else
            assigned_LUT{i} = GEOM.LUT_EPI; % Fallback
        end
    end
    
    % 4. Return structured results containing both times and assigned LUTs
    result = struct();
    result.init_dep = init_dep(:);
    result.init_rep = init_rep(:);
    result.assigned_LUT = assigned_LUT;
    
    % Minimal logging
    fprintf('FIXEDREPOLARIZATIONSCAN | T-wave peak (fixed REP): %.1f ms | DEP min/max: %.1f / %.1f ms\n', ...
            t_time, min(init_dep), max(init_dep));
end