function [parms, iter] = quamin(x, y, pinit, rfunc_gen)
    % Levenberg–Marquardt (jak quamin_pvd) dla dowolnego modelu owijanego
    % przez rfunc_gen (czyli używa Twojego generatora).

    arguments
        x 
        y 
        pinit 
        rfunc_gen function_handle = @rfunc
    end

    % START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [yest, G] = rfunc_gen(x, pinit, 1, []);

    res = y - yest;
    normresiter = norm(res);
    iter = 0; lambda = 0;
    ik = 1;
    RESNOW(ik,:) = [iter lambda normresiter];
    testnorm = normresiter;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loop1 = 1;
    parms = pinit(:).';          % zadbajmy o wektor wierszowy
    npar = length(pinit);
    if ~exist('noneg','var') == 1
        noneg = 0;
    end
    lambda = 0.001;

    % start iteration; outer loop
    while loop1 >= 1
        iter = iter+1;
        if iter > 1000
            break;
        end

        % compute GTG and G'*res
        GTG = G'*G;
        gtres = G'*res;
        break1 = 0;

        % compute new estimate
        while break1 == 0
            lambda = 2*lambda;
            if lambda > 1.e+10
                break1 = 1;
                break;
            end

            % Levenberg–Marquardt step
            MAT = GTG + lambda^2*eye(npar);
            delp = pinv(MAT) * gtres;

            if noneg == 1
                testp = max(parms + delp', 0);    % proste ograniczenie na >= 0
            else
                testp = parms + delp';
            end

            % test nowego wektora parametrów
            yest = rfunc_gen(x, testp, 0, []);
            res = y - yest;
            testnorm = norm(res);

            if testnorm < normresiter
               break;
            end
        end

        relres = abs(testnorm - normresiter) / normresiter;
        if relres < 1.e-5
            break1 = 1;
        end
        if break1 == 1
            break
        end

        parms = testp;
        normresiter = testnorm;
        lambda = lambda / 4;

        [~, G] = rfunc_gen(x, parms, 1, []);
    end
end
