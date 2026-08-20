function result = simplyfocusscan(GEOM, velocity)
    % SIMPLYFOCUSSCAN - Simple initialization of depolarization and repolarization times
    % based on RMS(BSM) peaks analysis and distance from the starting point.
    
    % Default wave propagation velocity (m/s)
    if nargin < 2 || isempty(velocity)
        velocity = 1; 
    end

    % 1. Calculate the RMS signal from all leads (BSM)
    rms_bsm = rms(GEOM.BSM, 1);
    
    % 2. Find peaks for the R-wave and T-wave
    % Find the global maximum (R-peak)
    [~, r_time] = max(rms_bsm);
    
    % Find the second peak (T-wave) in the time window after the R-wave
    % Starting the search 100 ms after the R-peak to skip the valley
    search_start = min(length(rms_bsm), r_time + 100);
    [~, t_time_offset] = max(rms_bsm(search_start:end));
    t_time = search_start - 1 + t_time_offset;
    
    % 3. Determine the single starting point (focus)
    % Project the potentials from r_time back onto the heart surface 
    % using the transpose of the AMA transfer matrix
    pot_heart = GEOM.AMA' * GEOM.BSM(:, r_time);
    
    % Select the node with the maximum absolute value as the starting point
    [~, focus_node] = max(pot_heart);
    
    % 4. Calculate depolarization times based on distance
    % Velocity: 1 m/s = 1 mm/s = 1 mm/ms
    vel_mm_ms = velocity; 
    
    % Euclidean distance from the starting point for each node
    diff_coords = GEOM.VER - GEOM.VER(focus_node, :);
    distances = sqrt(sum(diff_coords.^2, 2));
    
    % Depolarization time (delay relative to the starting point)
    dep = distances / vel_mm_ms;
    
    % 5. Calculate repolarization times
    rt = t_time - r_time;
    alpha = 0.8;
    rep = dep + (rt * alpha);
    
    % 6. Return results to the struct
    result = struct();
    result.dep = dep;
    result.rep = rep;
    result.r_time = r_time;
    result.t_time = t_time;
    result.focus_node = focus_node;
    
    % Minimal logging: min and max values for dep and rep
    fprintf('simplyfocusscan | DEP min/max: %4.1f / %4.1f ms | REP min/max: %4.1f / %4.1f ms | RT: %d ms\n', ...
            min(dep), max(dep), min(rep), max(rep), rt);
end