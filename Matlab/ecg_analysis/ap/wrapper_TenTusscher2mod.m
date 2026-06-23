function [t, V] = wrapper_TenTusscher2mod(HT, STOPTIME, CT, phases_mod, HR, ISO)
% WRAPPER_TENTUSSCHER2MOD Wrapper for the Ten Tusscher 2 model
% Generates several beats for a given HR to simulate the cell 
% adaptation phenomenon (restitution), and returns ONLY the LAST potential.

    % --- Handling optional arguments ---
    if nargin < 5 || isempty(HR)
        HR = 60;
    end

    if nargin < 6 || isempty(ISO)
        ISO = 0;
    end

    % --- Pacing parameters ---
    % Number of beats
    N_beats = 5; 
    
    % Heart cycle length in milliseconds
    BCL = 60000 / HR; 
    
    % Impulse amplitude and duration
    Stim_I = ones(1, N_beats) * -52;   
    Stim_T = ones(1, N_beats) * 1;     
    
    % Intervals: first impulse after 1 ms, subsequent ones every given BCL
    Stim_Int = [1, ones(1, N_beats - 1) * BCL]; 
    
    % Total simulation time:
    % Must accommodate N-1 BCL intervals and the time (STOPTIME) for the last beat.
    % +2 ms safety margin for initial offsets.
    TOTAL_STOPTIME = (N_beats - 1) * BCL + STOPTIME + 2;

    % --- Main model function call ---
    [t_full, V_full] = TenTusscher2mod(HT, TOTAL_STOPTIME, Stim_I, Stim_T, ISO, Stim_Int, CT, phases_mod);

    % --- Extraction of the last beat ---
    % Base time when the last beat started
    last_beat_start = (N_beats - 1) * BCL;
    
    % Find indices corresponding to the last run
    idx = find(t_full >= last_beat_start);
    
    % Reset the time vector so the last beat starts from ~0 ms again
    t_last = t_full(idx) - last_beat_start;
    V_last = V_full(idx);
    
    % Ensuring that the returned vector aligns EXACTLY with the STOPTIME window
    idx_stop = find(t_last <= STOPTIME);
    
    t = t_last(idx_stop);
    V = V_last(idx_stop);

end