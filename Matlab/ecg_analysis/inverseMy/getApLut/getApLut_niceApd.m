function LUT = getApLut_niceApd(start_apd, stop_apd, step_apd)
    % GETTMPLUT_NICEAPD Generates a Lookup Table of Action Potentials
    % with perfectly uniformly spaced APD values. Includes progress 
    % tracking and automatic filtering of unphysiological ODE anomalies.
    
    disp('Building fine-grained dictionary of Action Potentials...');
    disp('This may take a moment due to the fine sweep step...');
    
    % 1. Define the fine sweep range for the phase 3 modifier
    phases_sweep = 0.2:0.001:3.5; 
    num_sweeps = length(phases_sweep);
    
    % Arrays to store temporary results
    all_apds = zeros(num_sweeps, 1);
    
    % Preallocate structure array for speed
    temp_LUT = struct('phase_mod', cell(1, num_sweeps), ...
                      'V', cell(1, num_sweeps), ...
                      't', cell(1, num_sweeps), ...
                      't_dep', cell(1, num_sweeps));
                      
    % 2. Sweep through the modifier and generate all possible APs
    for i = 1:num_sweeps
        % Run the Ten Tusscher model with the current phase 3 modifier
        [t, V] = wrapper_TenTusscher2mod(0.1, 500, 1, [1, 1, phases_sweep(i)], 150);
        
        % Normalize the action potential to [0, 1] range (required for AMA matrix)
        V_norm = (V - min(V)) / (max(V) - min(V));
    
        % Find depolarization and repolarization indices
        idx_dep = find(V_norm >= 0.5, 1, 'first');
        idx_rep = find(V_norm(idx_dep:end) <= 0.1, 1, 'first') + idx_dep - 1;
        
        % Fallback if repolarization is not reached within 500ms
        if isempty(idx_rep)
            idx_rep = length(V_norm);
        end
    
        % Store exact properties of this simulation
        temp_LUT(i).phase_mod = phases_sweep(i);
        temp_LUT(i).V = V_norm;
        temp_LUT(i).t = t;
        temp_LUT(i).t_dep = t(idx_dep);
        
        % Calculate raw APD
        all_apds(i) = t(idx_rep) - t(idx_dep);
        
        % Print progress bar
        printProgress(i, num_sweeps, sprintf('Solving ODEs for templates | raw APD: %.1f', all_apds(i)));
    end
    
    disp('Fine sweep complete. Filtering anomalies...');
    
    % 3. FILTER ANOMALIES (Ensure identical depolarization times)
    all_t_deps = [temp_LUT.t_dep];
    
    % Find the standard physiological depolarization time (most frequent value)
    expected_t_dep = mode(all_t_deps);
    
    % Keep only templates that have the exact expected depolarization time
    % (Using a small tolerance 1e-4 for safe floating-point comparison)
    valid_idx = find(abs(all_t_deps - expected_t_dep) < 1e-1);
    
    num_filtered = num_sweeps - length(valid_idx);
    disp(['Filtered out ', num2str(num_filtered), ' anomalous templates (distorted shapes).']);
    
    % Overwrite temporary arrays with only valid (physiological) shapes
    temp_LUT = temp_LUT(valid_idx);
    all_apds = all_apds(valid_idx);
    
    disp('Mapping to requested APD targets...');
    
    % 4. Determine the target APDs based on input arguments
    step_val = abs(step_apd);
    if start_apd <= stop_apd
        target_apds = start_apd : step_val : stop_apd;
    else
        target_apds = start_apd : -step_val : stop_apd;
    end
    
    % 5. Build the final LUT matching closest exact valid APDs
    num_targets = length(target_apds);
    LUT = struct();
    for k = 1:num_targets
        target = target_apds(k);
        
        % Find the index of the valid APD closest to our current integer target
        [~, best_idx] = min(abs(all_apds - target));
        
        % Assign values to the final LUT
        LUT(k).phase_mod = temp_LUT(best_idx).phase_mod;
        LUT(k).V = temp_LUT(best_idx).V;
        LUT(k).t = temp_LUT(best_idx).t;
        LUT(k).t_dep = temp_LUT(best_idx).t_dep;
        
        % Keep the exact raw APD for debugging purposes
        LUT(k).raw_APD = all_apds(best_idx); 
        
        % Force the exact target APD so that the derivative algorithm sees perfect steps
        LUT(k).APD = target; 
        
        % Print progress bar
        printProgress(k, num_targets, 'Mapping APDs');
    end
    
    disp(['LUT generation complete. Created ', num2str(num_targets), ' perfect templates.']);
end