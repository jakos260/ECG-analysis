% function [score,notchopt,TST] = optimizeNotch(INV,notchopt,depopt,repopt,amplopt)
%
% Depending on the input, a certain ...
%
% INPUT:
% INV       = Inverse structure (inversefunc.m), which includes:
%               INV.T:          [matrix of nr of vertices on ventricles mesh x duration BSM]
%               INV.AMA:        The A-matrix
%               INV.ATA:        Transpose-A-matrix * A-matrix
%               INV.lpass:      Number of samples used in lowpassma filtering
%               INV.SPECS:      SPECS in INV-structure, comming from GEOM.SPECS [inversefunc.m]. 
%                               Contains inforamtion about slopes in RMS of signal. Needed to construct TMP's.
%               INV.mode:       The mode for getSfunc.m and getSmode.m [eg. depolarization, repolarization]
%               INV.PHIREF:     Measured BSM signals
%               INV.normphi:    Frobenius norm of measured BSM
%               INV.REGOP:      If INV.useSurfLapl == 1; REGOPREP = REGOP = normal surface laplacian
%               INV.REGOPREP:   ''
%               INV.ROTRO:      Transpose-surface laplacian-matrix * surface laplacian-matrix
%               INV.ROTROREP:   --
%               INV.lambdamax:  Criteria for ending iterations
%               INV.stopcrit:   Criteria for ending iterations
%
% notchopt  = structure containing notchopt.pot; zeros, unless INV.useNotch = 1 in inversefunc.m
% amplopt   = structure containing amplopt.pot; Amplitudes, defined in inversefunc.m as [OPT.AMP.pot = INV.amplitude]
%
% depopt    = structure with [.tims] and [.type]. These can be depolarization times if .type is 'dep' 
%             or repolarization times if .type is not 'dep'
% repopt    = a structure with [.tims]; reps if opt.tims are deps and vice versa
%
% OUTPUT:
% score     = [0] if new solution is not better than previous one, [1] if new solution is an improvement.
% notchopt  = same as input structure notchopt, but with new .pot [when score = 1]
% TST       = TST-structure with [PHIA, RES, rd, wrd, reg and tres].
%             These are the values that indicate how good your simulations fit the data in INV.PHIREF.
%
% Version 1: 01-apr-2015

function [score,notchopt,TST] = optimizeNotch(INV,notchopt,depopt,repopt,amplopt)

% prepare compution Marquardt step aimed at improving tims;
% compute Sprime and Sprime*Sprime' based on previous iteration
TST         = gettres_pot(INV,notchopt,depopt,repopt,amplopt);                                                      % Determine how good the initial input fits the data
starttres   = TST.tres;                                                                                             % Set start value of [RD plus regularization]
score       = 0;                                                                                                    % Set score initially to zero. Indicates if new .tims cause improvement of results.
lamb        = notchopt.lambopt;                                                                                     % Set lamb to initial lambdaopt value.
dif_pot     = 0.01;                                                                                                 % Change in .pot to produce Splus and Smin.
ampmod      = 4;                                                                                                    % mode for getSfunc.m

Splus       = getSfunc(INV.T, depopt.tims, repopt.tims, INV.SPECS, notchopt.pot+dif_pot, amplopt.pot, ampmod, INV); % Produce TMP's for .pot + dif_pot
Smin        = getSfunc(INV.T, depopt.tims, repopt.tims, INV.SPECS, notchopt.pot-dif_pot, amplopt.pot, ampmod, INV); % Produce TMP's for .pot - dif_pot
Sprime      = (Splus(:,notchopt.usetimes)-Smin(:,notchopt.usetimes)) / (dif_tims*2);
SST         = Sprime*Sprime';

% compute GTGM, with mu setting the weight of the regularization (smoothing based on the surface Laplacian)
M           = INV.ATA.*SST;
GTGM        = M'*M + + notchopt.mu^2*INV.ROTRO;                                                                     % why the + + and not just +??
gtres       = sum(INV.AMA.*(TST.RES*Sprime'))';                                                                     % Compute gtres
righths     = gtres - notchopt.mu^2*INV.ROTRO*notchopt.pot;

%% Compute next estimate; test residual;
% increase lambda it until starttres < TST.tres

while TST.tres >= starttres,
    
    lamb = lamb*2;
    
    if lamb > INV.lambdamax, break; end
    
    GTGML               = GTGM + lamb^2*eye(size(INV.AMA,2));
    delnotch            = inv(GTGML)*righths;
    newpot              = notchopt.pot + delnotch;                                                                	% Weight function to limit plateau phase apmlitue to at least 85%    
    newpot(newpot<0)    = 0;
    newpot(newpot>1)    = 1;
    testopt             = notchopt;
    testopt.pot         = newpot;
    TST                 = gettres_pot(INV,testopt,depopt,repopt,amplopt);                                           % Determine how good the new input [testopt] fits the data
    
    if TST.tres < starttres,                                                                                        % Check if new residual is smaller than in beginning
        if (starttres-TST.tres)/starttres >= INV.stopcrit,                                                          % Check if the improvement is big enough by comparing it to a stop criterion
            notchopt.pot    = newpot;                                                                              	% Change notchopt.pot to newpot, because of improvement
            score           = 1;                                                                                    % Change score to 1, because of improvement
            disp('optimize notch')
        end
        
        break; 
    end
end

% Set opt.lambopt to lamb devided by four [see AvO linsyst.pdf]:
notchopt.lambopt    = lamb/4;
TST                 = gettres_pot(INV,notchopt,depopt,repopt,amplopt);                                              % Calculate final result TST structure for output