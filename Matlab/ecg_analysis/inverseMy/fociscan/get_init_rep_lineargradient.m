function result = get_init_rep_lineargradient(GEOM, direction_vector, velocity)
    % REPOLARIZATIONGRADIENTSCAN - Initializes repolarization times and 
    % assigns specific TMP LUTs (EPI/ENDO) based on the endoVER flags.

    % Default wave propagation velocity: 1 m/s = 1 mm/ms
    if nargin < 3 || isempty(velocity)
        velocity = 1.0; 
    end
    
    % Check if direction vector is provided and valid
    if nargin < 2 || isempty(direction_vector) || numel(direction_vector) ~= 3
        error('A valid 1x3 or 3x1 direction_vector must be provided.');
    end
    
    % Ensure direction_vector is a column vector and normalize it
    u = direction_vector(:);
    u = u / norm(u);
    
    % 1. Calculate the spatial projection of each node onto the vector
    projections = GEOM.VER * u;
    distances = projections - min(projections);
    
    % 2. Calculate initial dep for rep calculations (1.0 m/s = 1.0 mm/ms)
    dep = distances / velocity; 
    
    % 3. Calculate repolarization times using BSM RMS peaks
    rms_bsm = rms(GEOM.BSM, 1);
    [~, r_time] = max(rms_bsm);
    
    search_start = min(length(rms_bsm), r_time + 100);
    [~, t_time_offset] = max(rms_bsm(search_start:end));
    t_time = search_start - 1 + t_time_offset;
    
    rt = t_time - r_time;
    alpha = 0.8;
    init_rep = dep + (rt * alpha);
    
    % 4. Assign potentials (LUTs) based on GEOM.endoVER
    num_nodes = length(GEOM.VER);
    assigned_LUT = cell(num_nodes, 1);
    
    for i = 1:num_nodes
        if GEOM.endoVER(i) == 0
            assigned_LUT{i} = GEOM.LUT_EPI;
        elseif GEOM.endoVER(i) == 1 || GEOM.endoVER(i) == 2
            assigned_LUT{i} = GEOM.LUT_ENDO;
        else
            assigned_LUT{i} = GEOM.LUT_EPI; % Fallback
        end
    end
    
    % 5. Return structured results
    result = struct();
    result.init_rep = init_rep;
    result.assigned_LUT = assigned_LUT;
    
    % Minimal logging
    fprintf('REPOLARIZATIONGRADIENTSCAN | Vel: %.2f m/s | REP min/max: %4.1f / %4.1f ms\n', ...
            velocity, min(init_rep), max(init_rep));
end