function [score,opt,TST] = my_optimizeDepRep(INV,opt,keepopt,notchopt,amplopt)
% MY_OPTIMIZEDEPREP - Corrected optimization function for dep/rep times.

TST         = gettres_v_from_ApLut(INV,opt,keepopt,notchopt,amplopt);
starttres   = TST.tres;
score       = 0;
lamb        = opt.lambopt;
dif_tims    = 1;

% 1. Determine times for derivative calculation
if strcmp(opt.type,'dep')
    tims_1 = opt.tims + dif_tims ;  tims_2 = opt.tims - dif_tims ;  
    tims_3 = keepopt.tims;          tims_4 = keepopt.tims;
else
    tims_1 = keepopt.tims;          tims_2 = keepopt.tims;          
    tims_3 = opt.tims + dif_tims ;  tims_4 = opt.tims - dif_tims ;
end

% 2. Calculate temporal derivative (Jacobian approximation)
Splus = my_regional_Smode(INV, tims_1, tims_3);
Smin  = my_regional_Smode(INV, tims_2, tims_4);

startopt = opt;
if (strcmp(opt.type,'dep') + strcmp(opt.type,'rep')) == 0
    startopt.tims = startopt.tims - keepopt.tims; 
end

Sprime = (Splus - Smin) / (dif_tims * 2);
SST    = Sprime * Sprime';

% 3. Compute regularization matrix
if strcmp(startopt.type,'dep')
    ROTRO_IN = INV.ROTRO; 
else
    ROTRO_IN = INV.ROTROREP; 
end
muREG = bsxfun(@times, opt.mu^2, ROTRO_IN);

% CRITICAL MATH FIX: Correct Levenberg-Marquardt Hessian and Gradient
H = INV.ATA .* SST; 
GTGM = H + muREG; 
gtres = sum(INV.AMA' .* (Sprime * TST.RES'), 2); 
righths = gtres - muREG * startopt.tims; 

testopt = startopt;
stuck = 0;
maxstuck = 10;

% 4. Iterate to find optimal LM damping parameter
while 1
    lamb = lamb * 2;
    
    if lamb > INV.lambdamax
        disp(['Optimization break: lambda reached maximum limit, lambda=', num2str(lamb)]); 
        break; 
    end
    
    % Damping factor for Levenberg-Marquardt
    GTGM_L = GTGM + lamb^2 * eye(size(INV.AMA, 2));
    
    % Solve for step direction
    deltau = GTGM_L \ righths;
    
    newtime = startopt.tims + deltau;
    % =========================================================================
    % PHYSIOLOGICAL CLAMPING (BOUNDS)
    % Forcing the optimizer to stay within realistic biology and LUT bounds
    % =========================================================================
    if strcmp(startopt.type, 'dep')
        % Depolarization cannot be negative (before QRS) and should not exceed e.g., 250 ms
        newtime = max(0, min(newtime, 150));
    elseif strcmp(startopt.type, 'apd')
        % APD MUST stay strictly within the bounds of your LUT!
        % (Assuming you generated LUT with start=200, stop=350)
        newtime = max(200, min(newtime, 460));
    end
    % =========================================================================

    testopt.tims = newtime;
    
    if strcmp(startopt.type,'apd')
        testopt.tims = testopt.tims + keepopt.tims; 
    end
    
    TST = gettres_v_from_ApLut(INV, testopt, keepopt, notchopt, amplopt);
    
    if TST.tres < starttres
        if (starttres - TST.tres) / starttres >= INV.stopcrit
            opt.tims = testopt.tims;
            score = 1;
            break;
        else
            stuck = stuck + 1;
            if stuck >= maxstuck
                disp(['Optimization break: stuck parameter reached ', num2str(stuck)]); 
                break; 
            end
        end
    else
        stuck = 0;
    end
end

opt.lambopt = lamb / 4;

TST = gettres_v_from_ApLut(INV, opt, keepopt, notchopt, amplopt);

end