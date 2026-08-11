function LUT = getApLut_HRmod(start_hr, stop_hr, step_hr)
    % GETTMPLUT_HRMOD Generates a Lookup Table of Action Potentials
    % with APD modulated by an increasing Heart Rate (HR).
    % Utilizes wrapper_TenTusscher2mod to handle the pacing train 
    % and extraction of the final steady-state beat automatically.
    
    disp('Building Action Potential Lookup Table based on HR sweep using wrapper...');
    
    % 1. Define the HR sweep range
    step_val = abs(step_hr);
    if start_hr <= stop_hr
        target_hrs = start_hr : step_val : stop_hr;
    else
        target_hrs = start_hr : -step_val : stop_hr;
    end
    
    num_targets = length(target_hrs);
    
    % Preallocate structure array for speed and memory optimization
    temp_LUT = struct('HR', cell(1, num_targets), ...
                      'CL', cell(1, num_targets), ...
                      'V', cell(1, num_targets), ...
                      't', cell(1, num_targets), ...
                      'APD', cell(1, num_targets), ...
                      't_dep', cell(1, num_targets));
                      
    valid_idx = false(1, num_targets);
    
    % 2. Sweep through requested Heart Rates
    for i = 1:num_targets
        current_hr = target_hrs(i);
        
        % Convert HR (beats per minute) to Cycle Length (CL) in milliseconds
        CL = 60000 / current_hr;
        
        % Wrapper configuration
        HT = 0.1;               % Integration step (ms)
        CT = 1;                 % Cell type: 1 = Epicardial
        phases_mod = [1, 1, 1]; % Neutral phase modifiers (no artificial APD scaling)
        ISO = 0;                % No sympathetic activation
        
        % STOPTIME handles the length of the extracted final beat.
        % 500 ms is standard, but we cap it at 95% of CL to prevent 
        % overshooting the physiological window at very high heart rates.
        STOPTIME = min(500, CL * 0.95);
        
        % Execute the wrapper simulation
        [t, V] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, phases_mod, current_hr, ISO);
        
        % Guard against solver crashes or truncated time arrays
        if isempty(t) || length(t) < 50
            continue; 
        end
        
        % Normalize the action potential to [0, 1] range (required for AMA matrix)
        V_norm = (V - min(V)) / (max(V) - min(V));
        
        % Find depolarization and repolarization indices
        idx_dep = find(V_norm >= 0.5, 1, 'first');
        if isempty(idx_dep)
            continue; % Skip non-depolarizing anomalies
        end
        
        idx_rep = find(V_norm(idx_dep:end) <= 0.1, 1, 'first') + idx_dep - 1;
        
        % Fallback if repolarization is not reached within the STOPTIME window
        if isempty(idx_rep)
            idx_rep = length(V_norm);
        end
        
        % Store exact properties of this simulation
        temp_LUT(i).HR = current_hr;
        temp_LUT(i).CL = CL;
        temp_LUT(i).V = V_norm;
        temp_LUT(i).t = t;
        temp_LUT(i).t_dep = t(idx_dep);
        temp_LUT(i).APD = t(idx_rep) - t(idx_dep);
        
        valid_idx(i) = true;
        
        % Print progress bar/logs
        if exist('printProgress', 'file') == 2
            printProgress(i, num_targets, sprintf('Mapping restitution via wrapper (HR = %.1f)', current_hr));
        else
            fprintf('Solved ODEs via wrapper for HR = %5.1f bpm (%d/%d)\n', current_hr, i, num_targets);
        end
    end
    
    disp('Sweep complete. Filtering anomalies...');
    
    % 3. Overwrite temporary arrays with only valid (physiological) shapes
    LUT = temp_LUT(valid_idx);
    
    num_filtered = num_targets - length(LUT);
    if num_filtered > 0
        disp(['Filtered out ', num2str(num_filtered), ' anomalous templates (distorted shapes).']);
    end
    
    disp(['LUT generation complete. Created ', num2str(length(LUT)), ' perfect templates.']);
end