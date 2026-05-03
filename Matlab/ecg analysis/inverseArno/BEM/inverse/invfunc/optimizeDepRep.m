% function [score,opt,TST] = optimizeDepRep(INV,opt,keepopt,notchopt,amplopt)
%
% Function to find new depolarization or repolarization times that better match with the measured signal.
%
% This function uses getSfunc.m to construct the TMP's.
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
% opt       = structure with [.tims], [.type], [.mu] and [.lambopt]. These can be depolarization times 
%             if .type is 'dep' or repolarization times if .type is not 'dep'.
% keepopt   = a structure with [.tims]; reps if opt.tims are deps and vice versa
% notchopt  = structure containing notchopt.pot; zeros, unless INV.useNotch = 1 in inversefunc.m
% amplopt   = structure containing amplopt.pot; Amplitudes, defined in inversefunc.m as [OPT.AMP.pot = INV.amplitude]
%
% OUTPUT:
% score     = [0] if new solution is not better than previous one, [1] if new solution is an improvement.
% opt       = same as input structure opt, but with new .tims [when score = 1]
% TST       = TST-structure with [PHIA, RES, rd, wrd, reg and tres].
%             These are the values that indicate how good your [dep/rep] times fit the data in INV.PHIREF.
%
% Version 1: 01-apr-2015
% New version JvdW:
% replaced getSfunc with getSmode so that new version of getSmode is used
% (see documentation on getSmode)
function [score,opt,TST] = optimizeDepRep(INV,opt,keepopt,notchopt,amplopt)

% Prepare compution Marquardt step aimed at improving tims;
% Compute Sprime and Sprime*Sprime' based on previous iteration

TST         = gettres_v(INV,opt,keepopt,notchopt,amplopt);                                                      % Determine how good the initial input fits the data
starttres   = TST.tres;                                                                                         % Set start value of [RD plus regularization]
score       = 0;                                                                                                % Set score initially to zero. Indicates if new .tims cause improvement of results.
lamb        = opt.lambopt;                                                                                      % Set lamb to initial lambdaopt value.
dif_tims    = 1;                                                                                                % Change in .tims to produce Splus and Smin.

%% Determine S for given .tims:

% Check if type is depolarization or something else ('rep' or 'apd'):
if strcmp(opt.type,'dep')
    tims_1 = opt.tims + dif_tims ;  tims_2 = opt.tims - dif_tims ;  tims_3 = keepopt.tims;          tims_4 = keepopt.tims;
else
    tims_1 = keepopt.tims;          tims_2 = keepopt.tims;          tims_3 = opt.tims + dif_tims ;  tims_4 = opt.tims - dif_tims ;
 end

% Calculate Splus and Smin:
% Splus      	= getSfunc(INV.T, tims_1, tims_3, INV.SPECS, notchopt.pot, amplopt.pot, INV.mode, INV);             % Produce TMP's for .tims + dif_tims
% Smin      	= getSfunc(INV.T, tims_2, tims_4, INV.SPECS, notchopt.pot, amplopt.pot, INV.mode, INV);             % Produce TMP's for .tims - dif_tims
Splus       = getSmode(INV.T, tims_1, tims_3, INV.SPECS, INV.mode, 22);    
Smin        = getSmode(INV.T, tims_2, tims_4, INV.SPECS, INV.mode, 22);
startopt   	= opt;                                                                                              % File initial ops [input] as startopt.

if (strcmp(opt.type,'dep') + strcmp(opt.type,'rep')) == 0, startopt.tims = startopt.tims - keepopt.tims; end    % In other words, if opt.type='apd'(?)

Sprime      = (Splus-Smin)/(dif_tims*2);                                                                        % Calculate average TMP between Splus and Smin
SST         = Sprime*Sprime';                                                                                   % ?

% Compute GTGM, with mu setting the weight of the regularization (smoothing based on the surface Laplacian):
if strcmp(startopt.type,'dep'), ROTRO_IN = INV.ROTRO; else; ROTRO_IN = INV.ROTROREP; end
muREG       = bsxfun(@times,opt.mu^2,ROTRO_IN);                                                                 % ?
M           = INV.ATA.*SST;                                                                                     % ?
GTGM        = M'*M + muREG;                                                                                     % ?
gtres       = sum(INV.AMA'.*(Sprime * TST.RES'),2);                                                             % original PvD
righths     = M'*gtres - muREG*startopt.tims;                                                                   % PvD original
testopt     = startopt;                                                                                         % File initial ops [input] [== startopt] as first testopt.

%% Compute next estimate; test residual;
% Increase lambda it until starttres < TST.tres
stuck       = 0;                                                                                                % number of consecutive iterations with small improvements
maxstuck    = 10;                                                                                               % Maximum number of tries without improvement of TST.tres in loop

% Perform this loop till stuck:
while 1
    lamb    = lamb*2;                                                                                           % Double lamb
    
    if lamb > INV.lambdamax, disp('lambdamax break'); break; end                                                % Check if lambdamax is still smaller than lambdamax-threshold, if not break from loop!
    
    GTGM_L          = GTGM + lamb^2*eye(size(INV.AMA,2));                                                       % Orginal AvO: see AvO linsyst.pdf
    deltau          = GTGM_L\righths;                                                                           % INVERSE; total computing time 0.87 times preceding
    newtime         = startopt.tims + deltau;                                                                   % Construct new times
    testopt.tims    = newtime;                                                                                  % Change .tims into new times.
    
    if strcmp(startopt.type,'apd'), testopt.tims = testopt.tims + keepopt.tims; end                             % ?
    
    TST             = gettres_v(INV,testopt,keepopt,notchopt,amplopt);                                          % Determine how good the new input [testopt.tims] fits the data
    
    if TST.tres < starttres                                                                                    % Check if new residual is smaller than in beginning
        if (starttres-TST.tres)/starttres >= INV.stopcrit                                                      % Check if the improvement is big enough by comparing it to a stop criterion
            opt.tims    = testopt.tims;                                                                         % Change opt.tims to testopt.tims, because of improvement
            score       = 1;                                                                                    % Change score to 1, because of improvement
            break                                                                                               % Break from loop
        else
            stuck       = stuck + 1;                                                                            % Improvement was not big enough, but do not stop procedure. Increase stuck parameter with one.
            
            if stuck >= maxstuck, fprintf('break stuck = %d\n',stuck); break; end                               % Check if stuck parameter has reached the maximum. If so, break from loop.
        end
    else
        stuck = 0;                                                                                              % still optimizing lambda
    end
end

% Set opt.lambopt to lamb devided by four [see AvO linsyst.pdf]:
opt.lambopt = lamb/4;