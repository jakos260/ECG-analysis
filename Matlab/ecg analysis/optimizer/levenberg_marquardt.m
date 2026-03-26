function [theta_opt, final_res_norm, iter] = levenberg_marquardt(obj_func, theta0, LB, UB, max_iter, tol)
% LEVENBERG_MARQUARDT Solves non-linear least squares problems with box constraints.
%
% INPUTS:
%   obj_func - Function handle: [res, G] = obj_func(p, need_G)
%              Returns residual vector (res) and Jacobian matrix (G).
%   theta0   - Initial guess for parameters (1 x npar).
%   LB       - Lower bounds for parameters (1 x npar). Empty [] means no lower limit.
%   UB       - Upper bounds for parameters (1 x npar). Empty [] means no upper limit.
%   max_iter - Maximum number of iterations.
%   tol      - Relative tolerance for the residual norm to stop the loop.
%
% OUTPUTS:
%   theta_opt      - Optimized parameters.
%   final_res_norm - Final Euclidean norm of the residuals.
%   iter           - Number of iterations performed.

    % Ensure initial guess is a row vector
    parms = theta0(:).';
    npar = length(parms);
    
    % Handle empty bounds by setting them to infinity
    if isempty(LB)
        LB = -inf(1, npar);
    else
        LB = LB(:).';
    end
    
    if isempty(UB)
        UB = inf(1, npar);
    else
        UB = UB(:).';
    end
    
    % Initial clamp (just in case theta0 is outside the defined bounds)
    parms = min(max(parms, LB), UB);
    
    % Initial evaluation
    [res, G] = obj_func(parms, true);
    normresiter = norm(res);
    
    iter = 0;
    lambda = 0.001;
    
    while iter < max_iter
        iter = iter + 1;
        
        % Compute GTG and G'*res
        GTG = G' * G;
        gtres = G' * res;
        break1 = false;
        
        % Inner loop to find valid step (Levenberg-Marquardt damping)
        while true
            lambda = 2 * lambda;
            if lambda > 1e10
                break1 = true;
                break;
            end
            
            % Calculate step change
            MAT = GTG + (lambda^2) * eye(npar);
            delp = pinv(MAT) * gtres;
            
            % Apply box constraints (clamping parameters between LB and UB)
            testp = parms + delp';
            testp = min(max(testp, LB), UB);
            
            % Test the new parameters
            [res_test, ~] = obj_func(testp, false);
            testnorm = norm(res_test);
            
            % If error decreased, accept step and break inner loop
            if testnorm < normresiter
                break;
            end
        end
        
        % Check convergence criteria based on relative change
        relres = abs(testnorm - normresiter) / normresiter;
        if relres < tol || break1
            break;
        end
        
        % Update parameters and lambda for the next iteration
        parms = testp;
        normresiter = testnorm;
        lambda = lambda / 4;
        
        % Recompute Jacobian at the new point
        [res, G] = obj_func(parms, true);
    end
    
    % Prepare outputs
    theta_opt = parms;
    final_res_norm = normresiter;
end