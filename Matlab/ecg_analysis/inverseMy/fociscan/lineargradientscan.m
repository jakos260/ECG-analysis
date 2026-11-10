function result = lineargradientscan(GEOM, direction_vector, vel_or_range)
    % LINEARGRADIENTSCAN - Initializes depolarization and repolarization times
    % by simulating a planar wave propagating along a specified 3D vector.
    %
    % The 3rd parameter (vel_or_range) can be either:
    %   - a scalar: represents wave propagation velocity (m/s)
    %   - a 1x2 vector [min, max]: represents the exact time range for REP (ms)

    % Default wave propagation velocity: 1 m/s = 1 mm/ms
    if nargin < 3 || isempty(vel_or_range)
        vel_or_range = 1.0; 
    end
    
    % Check if direction vector is provided and valid
    if nargin < 2 || isempty(direction_vector) || numel(direction_vector) ~= 3
        error('A valid 1x3 or 3x1 direction_vector must be provided.');
    end
    
    % Determine operation mode based on the 3rd parameter size
    if numel(vel_or_range) == 1
        mode = 'velocity';
        velocity = vel_or_range;
    elseif numel(vel_or_range) == 2
        mode = 'range';
        rep_min = min(vel_or_range);
        rep_max = max(vel_or_range);
        velocity = 1.0; % Default velocity for DEP calculations in range mode
    else
        error('The 3rd parameter must be a scalar (velocity) or a 2-element vector [min, max] (range).');
    end
    
    % Ensure direction_vector is a column vector and normalize it
    u = direction_vector(:);
    u = u / norm(u);
    
    % 1. Calculate the spatial projection of each node onto the vector
    projections = GEOM.VER * u;
    
    % Shift projections so the first activated node is at distance 0
    distances = projections - min(projections);
    
    % 2. Calculate depolarization times based on linear distance
    dep = distances / velocity; 
    
    % 3. Calculate ECG peaks (common for both modes)
    rms_bsm = rms(GEOM.BSM, 1);
    [~, r_time] = max(rms_bsm);
    
    search_start = min(length(rms_bsm), r_time + 100);
    [~, t_time_offset] = max(rms_bsm(search_start:end));
    t_time = search_start - 1 + t_time_offset;
    
    % 4. Calculate repolarization times depending on the mode
    if strcmp(mode, 'velocity')
        % Original logic based on BSM RMS peaks
        rt = t_time - r_time;
        alpha = 0.8;
        rep = dep + (rt * alpha);
    else
        % Range logic - scale distances linearly to the [min, max] range
        max_dist = max(distances);
        if max_dist == 0
            rep = rep_min * ones(size(distances));
        else
            dist_norm = distances / max_dist;
            rep = rep_min + dist_norm * (rep_max - rep_min);
        end
    end
    
    % 5. Return results
    result = struct();
    result.dep = dep;
    result.rep = rep;
    result.r_time = r_time;
    result.t_time = t_time;
    result.vector = u';
    
    % Minimal logging
    if strcmp(mode, 'velocity')
        fprintf('LINEARGRADIENTSCAN | Vel: %.2f m/s | DEP min/max: %4.1f / %4.1f ms | REP min/max: %4.1f / %4.1f ms\n', ...
                velocity, min(dep), max(dep), min(rep), max(rep));
    else
        fprintf('LINEARGRADIENTSCAN | Range mode | DEP min/max: %4.1f / %4.1f ms | REP min/max: %4.1f / %4.1f ms\n', ...
                min(dep), max(dep), min(rep), max(rep));
    end
end