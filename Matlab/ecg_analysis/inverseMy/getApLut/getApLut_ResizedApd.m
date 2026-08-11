function LUT = getApLut_ResizedApd(start_apd, stop_apd, step_apd)
    % GETTMPLUT_RESIZEDAPD Generates an Action Potential dictionary (LUT)
    % by generating one base AP and temporally interpolating it
    % to the desired widths (50% - 50% of amplitude).
    
    % 1. Base simulation parameters
    dt = 0.1; % 1 ms = 10 samples (10 kHz sampling)
    sim_time = 1500; % Extended time (1500 ms) so that the signal is not clipped when stretching the AP
    
    % Call the base model for a rigid phase 3 modifier: [1, 1, 1]
    [t_base, V_base] = wrapper_TenTusscher2mod(dt, sim_time, 1, [1, 1, 1], 150);
    
    % Normalization to the [0, 1] range
    V_norm = (V_base - min(V_base)) / (max(V_base) - min(V_base));
    
    % 2. Determination of the base APD (from 50% upstroke to 50% repolarization)
    idx_dep = find(V_norm >= 0.5, 1, 'first');
    idx_rep = find(V_norm(idx_dep:end) <= 0.5, 1, 'first') + idx_dep - 1;
    
    if isempty(idx_rep)
        error('Error: The base AP did not repolarize to 50% within the specified simulation time.');
    end
    
    t_dep_base = t_base(idx_dep);
    base_apd = t_base(idx_rep) - t_dep_base;
    
    disp(['base AP generated, APD: ', num2str(base_apd), ' ms']);
    
    % 3. Preparation of the target APD vector
    step_val = abs(step_apd);
    if start_apd <= stop_apd
        target_apds = start_apd : step_val : stop_apd;
    else
        target_apds = start_apd : -step_val : stop_apd;
    end
    
    num_targets = length(target_apds);
    
    % Preallocate structure for performance
    LUT = repmat(struct('phase_mod', 1, 'V', [], 't', [], 't_dep', [], 'raw_APD', [], 'APD', []), 1, num_targets);
    
    % 4. Scaling and creating the final LUT
    for k = 1:num_targets
        target = target_apds(k);
        
        % Scaling factor (ratio of target APD to base APD)
        scale_factor = target / base_apd;
        
        % Time axis transformation. We anchor the t_dep_base point, 
        % so as not to shift the activation moment in time (the upstroke stays in place).
        t_stretched = t_dep_base + (t_base - t_dep_base) * scale_factor;
        
        % Resampling (interpolation) of the signal to the original, constant time grid (10 kHz)
        % We use the 'V_norm(end)' value for right-side extrapolation (to avoid NaNs)
        V_resized = interp1(t_stretched, V_norm, t_base, 'linear', V_norm(end));
        
        % Correction of numerical errors from interpolation
        V_resized(V_resized < 0) = 0;
        V_resized(V_resized > 1) = 1;
        
        % Finding new indices to verify the final width
        idx_dep_new = find(V_resized >= 0.5, 1, 'first');
        idx_rep_new = find(V_resized(idx_dep_new:end) <= 0.5, 1, 'first') + idx_dep_new - 1;
        
        if ~isempty(idx_rep_new)
            raw_apd = t_base(idx_rep_new) - t_base(idx_dep_new);
        else
            raw_apd = NaN;
        end
        
        % Save the results
        LUT(k).phase_mod = 1; 
        LUT(k).V = V_resized;
        LUT(k).t = t_base;
        LUT(k).t_dep = t_base(idx_dep_new);
        LUT(k).raw_APD = raw_apd; 
        LUT(k).APD = target;
        
        % Optional progress bar if you use a separate function in your workspace
        % printProgress(k, num_targets, sprintf('Creating APD: %.1f ms', target));
    end
    
    disp(['created ', num2str(num_targets), ' AP templates']);
end