function result = lineargradientscan(GEOM, direction_vector, velocity)
    % LINEARGRADIENTSCAN - Initializes depolarization and repolarization times
    % by simulating a planar wave propagating along a specified 3D vector.

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
    % GEOM.VER contains coordinates in millimeters (Nx3)
    projections = GEOM.VER * u;
    
    % Shift projections so the first activated node is at distance 0
    distances = projections - min(projections);
    
    % 2. Calculate depolarization times based on linear distance
    % Velocity conversion: 1.0 m/s = 1000 mm/s = 1.0 mm/ms
    vel_mm_ms = velocity; 
    
    dep = distances / vel_mm_ms;
    
    % 3. Calculate repolarization times using BSM RMS peaks (same as simplyfocusscan)
    rms_bsm = rms(GEOM.BSM, 1);
    [~, r_time] = max(rms_bsm);
    
    search_start = min(length(rms_bsm), r_time + 100);
    [~, t_time_offset] = max(rms_bsm(search_start:end));
    t_time = search_start - 1 + t_time_offset;
    
    rt = t_time - r_time;
    alpha = 0.8;
    rep = dep + (rt * alpha);
    
    % 4. Return results
    result = struct();
    result.dep = dep;
    result.rep = rep;
    result.r_time = r_time;
    result.t_time = t_time;
    result.vector = u';
    
    % Minimal logging
    fprintf('LINEARGRADIENTSCAN | Vel: %.2f m/s | DEP min/max: %4.1f / %4.1f ms | REP min/max: %4.1f / %4.1f ms\n', ...
            velocity, min(dep), max(dep), min(rep), max(rep));
end