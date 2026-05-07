% TST = gettres_v(INV, opt, keepopt, notchopt, amplopt)
%
% Determine how good the times [eg. opt.tims] fit your data [INV.PHIREF].
%
% INPUT:
% INV       = Inverse structure (inversefunc.m), which includes:
%               INV.T:          [matrix of nr of vertices on ventricles mesh x duration BSM]
%               INV.AMA:        The A-matrix
%               INV.lpass:      Number of samples used in lowpassma filtering
%               INV.SPECS:      SPECS in INV-structure, comming from GEOM.SPECS [inversefunc.m]. 
%                               Contains inforamtion about slopes in RMS of signal. Needed to construct TMP's.
%               INV.mode:       The mode for getSfunc.m and getSmode.m [eg. depolarization, repolarization]
%               INV.PHIREF:     Measured BSM signals
%               INV.normphi:    Frobenius norm of measured BSM
%               INV.REGOP:      If INV.useSurfLapl == 1; REGOPREP = REGOP = normal surface laplacian
%               INV.REGOPREP:   ''
%
% opt       = structure with [.tims] and [.type]. These can be
%             depolarization times if .type is 'dep' or repolarization times if .type is not 'dep'
% keepopt   = a structure with [.tims]; reps if opt.tims are deps and vice versa
% notchopt  = structure containing notchopt.pot; zeros, unless INV.useNotch = 1 in inversefunc.m
% amplopt   = structure containing amplopt.pot; Amplitudes, defined in inversefunc.m as [OPT.AMP.pot = INV.amplitude]
%
% OUTPUT:
% TST       = TST-structure with [PHIA, RES, rd, wrd, reg and tres].
%             These are the values that indicate how good your [dep/rep] times fit the data in INV.PHIREF.
%
% Version 1: 01-apr-2015
% New version JvdW:
% replaced getSfunc with getSmode so that new version of getSmode is used
% (see documentation on getSmode)

function TST = gettres_v(INV, opt, keepopt, notchopt, amplopt)

% Set weighted RD value parameter to 0.0010;
wrd_param   = 0.0010;

%% Determine S for given .tims:
% Check if type is depolarization or something else ('rep' or 'apd'):
if strcmp(opt.type,'dep'),
    tims_1 = opt.tims; tims_2 = keepopt.tims;
else
    tims_1 = keepopt.tims; tims_2 = opt.tims;
end

% TST.S = getSfunc(INV.T, tims_1, tims_2, INV.SPECS, notchopt.pot, amplopt.pot, INV.mode, INV);   %
TST.S = getSmode(INV.T, tims_1, tims_2, INV.SPECS, INV.mode); % 22: test new option with estimation based on tdom

%% Determine simulated BSM, residue, RD and weighted RD values:
TST.PHIA    = lowpassma(INV.AMA * TST.S,INV.lpass);                                             % Construct filterd simulated BSM with A-matrix (INV.AMA) and simulated cardiac signals (TST.S)
TST.RES     = INV.PHIREF-TST.PHIA(1:size(INV.PHIREF,1),1:size(INV.PHIREF,2));                   % Calculated residue --> difference INV.PHIREF and TST.PHIA
TST.rd      = norm(TST.RES,'fro')/INV.normphi;                                                  % Calculate RD value
TST.wrd     = sum(rms(TST.RES)./(wrd_param + rms(INV.PHIREF)));                                 % Calcualte weighted RD value

%% Determine regularization value:
if strcmp(opt.type,'dep'),
    TST.reg = norm(INV.REGOP*opt.tims);                                                         % Calculate regularization by multiplying REGOP with times
else
    TST.reg = norm(INV.REGOPREP*opt.tims);                                                      % Calculate regularization by multiplying REGOPREP with times
end

%% Determine value of regularization term plus RD value:
if INV.useWeighedRd                                                                            % Check if weighted RD is used:
    if strcmp(opt.type,'dep'),
        TST.tres = sqrt(TST.wrd^2+(TST.reg*opt.mu)^2);
    else
        TST.tres = sqrt(TST.wrd^2+(TST.reg*opt.mu)^2);
    end
else
    if strcmp(opt.type,'dep'),
        TST.tres = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
    else
        TST.tres = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);
    end
end