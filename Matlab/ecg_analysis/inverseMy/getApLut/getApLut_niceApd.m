function LUT = getApLut_niceApd(start_apd, stop_apd, step_apd, CT)
    disp('Building fine-grained dictionary of Action Potentials...');
    disp('This may take a moment due to the fine sweep step...');
    
    num_sweeps = 1500; % number of sweep samples
    phase2_sweep = linspace(0.5, 1.9, num_sweeps);
    phase3_sweep = linspace(3.0, 0.5, num_sweeps);
    
    all_apds = NaN(num_sweeps, 1);
    
    temp_LUT = struct('phase_mod', cell(1, num_sweeps), ...
                      'V', cell(1, num_sweeps), ...
                      't', cell(1, num_sweeps), ...
                      't_dep', cell(1, num_sweeps));
                      
    for i = 1:num_sweeps
        printProgress(i, num_sweeps, sprintf('Solving ODEs | raw APD: %.1f', all_apds(i)));

        phase2_mod = phase2_sweep(i);
        phase3_mod = phase3_sweep(i);
        
        sim_time = 800;
        
        [t, V] = wrapper_TenTusscher2mod(0.1, sim_time, CT, [1, phase2_mod, phase3_mod], 100);
        
        % ==========================================================
        % 1. Potential filtering before normalization
        % ==========================================================
        if V(1) > -75
            continue;
        end
        
        if (max(V) - min(V)) < 80
            continue;
        end
        
        V_norm = (V - min(V)) / (max(V) - min(V));
    
        idx_dep = find(V_norm >= 0.5, 1, 'first');
        if isempty(idx_dep); continue; end
        
        % ==========================================================
        % 2. Filtering "Loss of Dome" (for Epicardium)
        % ==========================================================
        if CT == 1
            % peak phase 0 occured aproximately 2-4 ms after idx_dep.
            % searching for notch between 5 ms and 40 ms after repolarization.
            % searching for plateau peak (dome) between 40 ms and 100 ms.
            idx_5ms   = idx_dep + 50;
            idx_40ms  = idx_dep + 400;
            idx_100ms = idx_dep + 1000;
            
            if idx_100ms <= length(V_norm)
                dome_peak = max(V_norm(idx_40ms:idx_100ms));
                notch_min = min(V_norm(idx_5ms:idx_40ms));
                
                if dome_peak < 0.75
                    continue;
                end
                
                if notch_min < 0.60
                    continue;
                end
            end
        end
        % ==========================================================
        
        idx_rep = find(V_norm(idx_dep:end) <= 0.1, 1, 'first');
        if isempty(idx_rep); continue; end
        idx_rep = idx_rep + idx_dep - 1;
        
        if max(V_norm(idx_rep:end)) > 0.15; continue; end
        
        % optional safety test
        idx_rest = find(V_norm(idx_rep:end) <= 0.05, 1, 'first');
        if isempty(idx_rest)
            continue; 
        end
   
        temp_LUT(i).phase_mod = [phase2_mod, phase3_mod];
        temp_LUT(i).V = V_norm;
        temp_LUT(i).t = t;
        temp_LUT(i).t_dep = t(idx_dep);
        
        all_apds(i) = t(idx_rep) - t(idx_dep);
    end
    
    disp('Fine sweep complete. Removing anomalies...');
    
    valid_sims = ~isnan(all_apds);
    temp_LUT = temp_LUT(valid_sims);
    all_apds = all_apds(valid_sims);
    
    if isempty(all_apds)
        error('Error: None of the simulations passed the filters! Refine the search mesh (phases_sweep) or loosen the filters for CT=%d.', CT);
    end
    
    all_t_deps = [temp_LUT.t_dep];
    expected_t_dep = mode(all_t_deps);
    
    valid_idx = find(abs(all_t_deps - expected_t_dep) < 5.0);
    
    num_filtered = length(all_apds) - length(valid_idx);
    disp(['Filtered out ', num2str(num_filtered), ' anomalous templates.']);
    
    temp_LUT = temp_LUT(valid_idx);
    all_apds = all_apds(valid_idx);
    
    disp('Mapping to requested APD targets...');
    
    step_val = abs(step_apd);
    if start_apd <= stop_apd
        target_apds = start_apd : step_val : stop_apd;
    else
        target_apds = start_apd : -step_val : stop_apd;
    end
    
    num_targets = length(target_apds);
    LUT = struct();
    for k = 1:num_targets
        target = target_apds(k);
        
        [min_diff, best_idx] = min(abs(all_apds - target));
        
        if min_diff > 3.0
            warning('Cel APD %.1f ms nieosiągalny. Najbliższa bezpieczna wartość: %.1f ms.', target, all_apds(best_idx));
        end
        
        LUT(k).phase_mod = temp_LUT(best_idx).phase_mod;
        LUT(k).V = temp_LUT(best_idx).V;
        LUT(k).t = temp_LUT(best_idx).t;
        LUT(k).t_dep = temp_LUT(best_idx).t_dep;
        LUT(k).raw_APD = all_apds(best_idx);
        LUT(k).APD = target;
        
        printProgress(k, num_targets, 'Mapping APDs');
    end
    
    disp(['LUT generation complete. Created ', num2str(num_targets), ' perfect templates.']);
end

%%
CT = 1;
LUT = getApLut_niceApd(150, 500, 1, CT);

filename = fullfile('inverseMy', 'getApLut', sprintf('ApLut_niceApd_CT%d.mat', CT));
save(filename, 'LUT');

clf
hold on
for i=1:size(LUT,2)
    plot(LUT(i).V);
end