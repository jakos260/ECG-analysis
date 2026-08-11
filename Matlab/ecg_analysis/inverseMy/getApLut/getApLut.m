function LUT = getApLut(start, stop, step)
    % phases_sweep = 0.8:0.02:1.2;
    phases_sweep = start:step:stop;
    LUT = struct();
    for i = 1:length(phases_sweep)
        [t, V] = wrapper_TenTusscher2mod(0.1, 500, 1, [1, 1, phases_sweep(i)], 60, 0);
        V_norm = (V - min(V)) / (max(V) - min(V));
    
        idx_dep = find(V_norm >= 0.5, 1, 'first');
        idx_rep = find(V_norm(idx_dep:end) <= 0.1, 1, 'first') + idx_dep - 1;
    
        LUT(i).phase_mod = phases_sweep(i);
        LUT(i).V = V_norm;
        LUT(i).APD = t(idx_rep) - t(idx_dep);
    end
end