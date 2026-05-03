% TST = gettres_pot(INV,opt,depopt,repopt,otheropt)
%
% Determine how good the ...
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
%               INV.REGOP:      If INV.useSurfLapl == 1; REGOPREP = REGOP = normal surface laplacian
%
% opt       =  amplopt  [in optimizeAmpl.m]
% otheropt  =  notchopt [in optimizeAmpl.m]
% notchopt  = structure containing notchopt.pot; zeros, unless INV.useNotch = 1 in inversefunc.m
% amplopt   = structure containing amplopt.pot; Amplitudes, defined in inversefunc.m as [OPT.AMP.pot = INV.amplitude]
%
% depopt    = structure with [.tims] and [.type]. These can be depolarization times if .type is 'dep' 
%             or repolarization times if .type is not 'dep'
% repopt    = a structure with [.tims]; reps if opt.tims are deps and vice versa
%
% OUTPUT:
% TST       = TST-structure with [PHIA, RES, rd, wrd, reg and tres].
%             These are the values that indicate how good your [dep/rep]
%             times fit the data in INV.PHIREF.
%
% Version 1: 01-apr-2015

function TST = gettres_pot(INV,opt,depopt,repopt,otheropt)

%% Determine S for given .tims:
% Check if type is notch or something else:
if strcmp(opt.type,'notch'),
    pot_1 = opt.pot; pot_2 = otheropt.pot;
else
    pot_1 = otheropt.pot; pot_2 = opt.pot;
end

S = getSfunc(INV.T, depopt.tims, repopt.tims, INV.SPECS, pot_1, pot_2, INV.mode, INV);

%% Determine simulated BSM, residue, RD and weighted RD values:
TST.PHIA    = lowpassma(INV.AMA*S(:,opt.usetimes),INV.lpass);               % Construct filterd simulated BSM with A-matrix (INV.AMA) and simulated cardiac signals (TST.S)
TST.RES     = INV.PHIREF(:,opt.usetimes) - TST.PHIA;                        % Calculated residue --> difference INV.PHIREF and TST.PHIA
TST.rd      = norm(TST.RES,'fro')/norm(INV.PHIREF(:,opt.usetimes),'fro');   % Calculate RD value
TST.reg     = norm(INV.REGOP*opt.pot);                                      % Calculate regularization by multiplying REGOP with times
TST.tres    = sqrt(TST.rd^2+(TST.reg*opt.mu)^2);                            % Calculate residue of RD plus regularization
TST.S       = S;                                                            % File TMP's for output.