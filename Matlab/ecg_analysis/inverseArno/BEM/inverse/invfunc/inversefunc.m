% function meas = inversefunc(varargin)
%
% Inverse determination of timing of Equivalent Double Layer source [EDL].
% The optimization is based on basic Marquardt type of solving the non-linear
% estimation problem. Alternate iterations between depolarization and
% repolarization is an option.
%
% OUTPUT:
% meas              = Structure with final depolarization times,
%                     repolarization times, RD value, REG value &
%                     correlation value.
%
% INPUT:
% varargin{1}       = GEOM (Structure created with InvInit.m.)
% varargin{2}       = Initial depolarization times {multifociscan.m}
% varargin{3}       = Initial repolarization times {multifociscan.m}
%
% OPTIONAL INPUT:
% 'mudep'           = regularization parameter for depolarization [default = 1.5e-4]
% 'murep',          = regularization parameter for repolarization [default = 1.5e-4]
% 'repopt',         = type of repolarization [default = 'apd']
% 'muampl',         = regularization parameter for depolarization [default = 1.5e-4]
% 'estimateampl',   = [1] use optimizeAmpl.m instead of optimizeDepRep.m [default = 0]
% 'maxiter',        = maximum number of iterations [default = 25]
% 'minrd',          = RD value to be reached [default = 0.15]
% 'casedir',        = directory for results [default = 'invresults\']
% 'logname',        = ... [default = '']
% 'leads',          = ... [default = 1:65]
% 'weighedrd',      = [1] Use weighed RD instead of normal RD [default = 0]
% 'reg',            = ... [default = [] ]
% 'blmode',         = [1] blcorrection qrsonset - T, [2] (Vstim-10) - T. [default = 1]
% 'lpass',          = moving avarage lowpass filtering. number of samples. [default = 5]
% 'mode',           = Standard mode: [1] depolarization [4] depolarization and repolarization [default = 4]
%
% INV.mode = 1, depolarization only
% INV.mode = 2, repolarization only
% INV.mode = 3, notch?
% INV.mode = 4, depolarization and repolarization
% INV.mode = 5, depolarization and repolarization & notch?
% INV.mode = 6, AP amplitude only?
%
% Peter van Dam; 2010 november.
% All rights reserved Peacs, Arnhem  the Netherlands
%
% Copy of inverse_po.m 9-jan-2015: version 1.0 AJ
%
% Version 1: 01-apr-2015

function meas = inversefunc(varargin)

%% inverse procedure parameters

% STANDARD PARAMETERS
INV.mudep           = 1.5e-4;           % regularization parameter for depolarization
INV.murep           = 1.5e-4;           % regularization parameter for repolarization
INV.muampl          = 1.5e-4;           % regularization parameter for ... ?
INV.repOpt          = 'rep';            % type of repolarization
INV.MINRD           = 0.15;             % minimum RD
INV.maxiter         = 25;               % maximum number of iterations
INV.lambopt         = 0.001;            % initial value of lambda in inverse procedure [e.g. optimizeDepRep.m]
INV.useAmpl         = 0;                % [1] use optimizeAmpl.m
INV.useNotch        = 0;                % [1] use optimizeNotch.m
INV.useSurfLapl     = 1;                % [1] Calculate normal surface laplacian with calcREGOP.m
INV.useWeighedRd    = 0;                % [1] Use weighted RD instead of normal RD
INV.casedir         = 'invresults\';    % directory for results
INV.doSave          = 0;                %
INV.subname         = '';               %
INV.leads           = 1:65;             %
INV.mode            = 4;                % Standard mode: [1] depolarization [4] depolarization and repolarization
INV.doAmplitude     = 0;                %
INV.reg             = [];               %
INV.blmode          = 1;                % [1] blcorrection qrsonset - T, [2] (Vstim-10) - T.
INV.lpass           = 5;                % moving avarage lowpass filtering. number of samples.
INV.lambdamax       = 500;              % criteria for ending iterations
INV.stopcrit        = 2e-4;             % criteria for ending iterations
versie              = 1;                %
min_Vstim           = 10;               % parameter to indicate the number of samples taken into account BEFORE stimulation marker
notch_amp_param_1   = 10;               % parameter used for amplitude and notch optimization
notch_amp_param_2   = 40;               % parameter used for amplitude and notch optimization
notch_amp_param_3   = 60;               % parameter used for amplitude and notch optimization
miniter             = 2;                % minimum number of iterations to be performed!

% ADJUST PARAMETERS WITH INPUT:
if length(varargin) < 3
    error('This routine needs at least three parameters');
else
    GEOM            = varargin{1};
    INV.amplitude   = ones(size(GEOM.VER,1),1);
    INV.initdep     = varargin{2};
    INV.initrep     = varargin{3};
    INV.leads       = 1:length(GEOM.TVER);
    usetimes        = size(GEOM.BSM,2) - GEOM.SPECS.onsetqrs;
    pp              = 4;
    while pp <= nargin
        if ischar(varargin{pp})
            key = lower(varargin{pp});
            switch key
                case 'mudep', INV.mudep              = varargin{pp+1}; pp = pp+2;
                case 'murep', INV.murep              = varargin{pp+1}; pp = pp+2;
                case 'repopt', INV.repOpt            = varargin{pp+1}; pp = pp+2;
                case 'muampl', INV.muampl            = varargin{pp+1}; pp = pp+2;
                case 'estimateampl', INV.useAmpl     = varargin{pp+1}; pp = pp+2;
                case 'maxiter', INV.maxiter          = varargin{pp+1}; pp = pp+2;
                case 'minrd', INV.MINRD              = varargin{pp+1}; pp = pp+2;
                case 'casedir', INV.casedir          = varargin{pp+1}; pp = pp+2;
                case 'logname', INV.subname          = varargin{pp+1}; pp = pp+2;
                case 'leads', INV.leads              = varargin{pp+1}; pp = pp+2;
                case 'weighedrd', INV.useWeighedRd   = varargin{pp+1}; pp = pp+2;
                case 'mode', INV.mode                = varargin{pp+1}; pp = pp+2;
                case 'reg', INV.reg                  = varargin{pp+1}; pp = pp+2;
                case 'blmode', INV.blmode            = varargin{pp+1}; pp = pp+2;
                case 'lpass', INV.lpass              = varargin{pp+1}; pp = pp+2;
                case 'doSave', INV.doSave            = varargin{pp+1}; pp = pp+2;
                case 'dosave', INV.doSave            = varargin{pp+1}; pp = pp+2;
                otherwise, error('unknown parameter');
            end
        end
    end
end

if ~strcmp(INV.repOpt,'rep'), INV.repOpt = 'apd'; end                                               % when INV.repOpt is not 'rep', change it to 'apd'
if ~exist(INV.casedir,'dir'), mkdir(INV.casedir); end                                               % check if output directory is present, else create directory

% Create names of files:
% subname = [ '_' num2str(size(GEOM.LAY,1)-1)];                                                       % make sub name with number of ECG electrodes in GEOM.LAY
runcase = fullfile(INV.casedir,[GEOM.subject '_' GEOM.beat]);                                   % name of file

%% PREPARE INPUT DATA
INV.VER     = GEOM.VER;
INV.ITRI    = GEOM.ITRI;
INV.ADJ     = GEOM.ADJ;
INV.DIST    = GEOM.DIST;
INV.neigh   = GEOM.neigh;

% BASELINE CORRECTION FOR DATA [BSM]:
if strcmp(GEOM.type,'_atria')
    
    if INV.blmode ~= 1, warning('MATLAB:UndefinedFunction', 'Only baselinemode 1 supported for atria'); end
    
    INV.BSM =  baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endtwave),INV.lpass));	% Baseline correction and low pass filtering for the atria
    
    if INV.mode == 1, usetimes_start = GEOM.SPECS.onsetqrs; else; usetimes_start = usetimes; end
    
else
    
    % [1] blcorrection qrsonset - T, [2] (Vstim - min_Vstim) - T:
    if INV.blmode == 1
        INV.BSM = baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),INV.lpass));
    else
        bsmt    = baselinecor(GEOM.BSM,GEOM.SPECS.time_Vstim - min_Vstim,GEOM.SPECS.endtwave);
        INV.BSM = lowpassma(bsmt(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),INV.lpass);
    end
    
    if INV.mode == 1, usetimes_start = GEOM.SPECS.qrsduration; else; usetimes_start = usetimes; end
    
end

INV.usetimes    = min(usetimes_start ,size(INV.BSM,2));                                     % Define usetimes --> end timestamp for analyses.
INV.PHIREF      = INV.BSM(:,1:INV.usetimes);                                                % Define reference signal from BSM.
INV.normphi     = norm(INV.PHIREF,'fro');                                                   % Frobenius norm of BSM (during QRS)
INV.T           = ones(size(GEOM.AMA,2),1)*(0:INV.usetimes-1);                              % matrix of nr of vertices on ventricles mesh x duration BSM (during QRS)

% Compute regularization operator, e.g.surface Laplacian:
if isempty(INV.reg)
    [INV.REGOP,INV.REGOPREP]    = calcREGOP(GEOM,INV.useSurfLapl);                          % IF INV.useSurfLapl == 1; REGOPREP = REGOP = normal surface laplacian
else
    INV.REGOP                   = INV.reg;
    INV.REGOPREP                = INV.REGOP;
end

% Precompute large matrixes:
INV.AMA         = GEOM.AMA;                                                                 % Transfer (A) Matrix
INV.ATA         = INV.AMA'*INV.AMA;                                                         % transpose-A-matrix * A-matrix
INV.ROTRO       = INV.REGOP'*INV.REGOP;                                                     % transpose-surface laplacian-matrix * surface laplacian-matrix
INV.ROTROREP    = INV.REGOPREP'*INV.REGOPREP;                                               % IF INV.useSurfLapl == 1; INV.ROTROREP = INV.ROTRO
INV.SPECS       = GEOM.SPECS;                                                               % File SPECS in INV-structure.

%% read initial estimates;
% Regularization parameters must be used to tune the final result: small values may produce spatially 'wild' solutions: use trial and error

% Depolarization
OPT.DEP.mu      = INV.mudep;                                                                % Regularization parameter for depolarization
OPT.DEP.lambopt = INV.lambopt;                                                              % ?
OPT.DEP.tims    = INV.initdep;                                                              % Initial depolarization times
OPT.DEP.type    = 'dep';                                                                    % Tag for type

% Repolarization
OPT.REP.mu      = INV.murep;                                                                % Regularization parameter for repolarization
OPT.REP.lambopt = INV.lambopt;                                                              % ?
OPT.REP.tims    = INV.initrep;                                                              % Initial repolarization times
OPT.REP.type    = INV.repOpt;                                                               % Tag for type 'rep' or 'apd'

% Jpoint:
OPT.NOT.pot     = zeros(size(GEOM.VER,1),1);
if INV.useNotch                                                                             % Future purpose; when also a notch is included in the source model:
    OPT.NOT.pot      = 1-sqrAmpl(GEOM,GEOM.SPECS.time_Jpoint+notch_amp_param_1:GEOM.SPECS.time_Jpoint+notch_amp_param_2); %
    OPT.AMP.usetimes = GEOM.SPECS.time_Jpoint:GEOM.SPECS.time_Jpoint+notch_amp_param_3;     %              %
    OPT.NOT.mu       = INV.muampl;                                                          % Regularization parameter for amplitude
    OPT.NOT.lambopt  = INV.lambopt;                                                         % ?
    OPT.NOT.type     = 'notch';                                                            	% Tag for type
end

% AP amplitude:
OPT.AMP.pot = INV.amplitude;
if INV.useAmpl
    OPT.AMP.pot      = sqrAmpl(GEOM,GEOM.SPECS.time_Jpoint+notch_amp_param_1:GEOM.SPECS.time_Jpoint+notch_amp_param_2); %
    OPT.AMP.usetimes = GEOM.SPECS.time_Jpoint:GEOM.SPECS.time_Jpoint+notch_amp_param_3;     %
    OPT.AMP.lambopt  = INV.lambopt;                                                         % ?
    OPT.AMP.mu       = INV.muampl;                                                          % Regularization parameter for amplitude
    OPT.AMP.type     ='amplitude';                                                        	% Tag for type
    INV.doAmplitude  = 1;                                                                   % Set to 1, for Amplitude analyses with optimizeAmpl.m
end

% Determine how good the initial estimate fits the data with gettres_v.m:
TST = gettres_v(INV, OPT.DEP, OPT.REP, OPT.NOT, OPT.AMP);

%% =======================================================================
% DOCUMENT CASE:
if INV.doSave
    logfile = [runcase '.log'];
    outfil  = [runcase '.src'];
    fh      = fopen(logfile,'w');
    fprintf(fh,'%s\n',['file: ' logfile]);
    fprintf(fh,'%s\n',['date: ' datestr(clock)]);
    fprintf(fh,'%s\n',['prog: invgeom (' num2str(versie) ')' ]);
    fprintf(fh,'%s\n',['initvals:   ' [runcase 'srcinit']]);
    fprintf(fh,'%s\n',['subject:   ' GEOM.subject]);
    fprintf(fh,'%s\n',['used anisotropyRatio:    ' GEOM.anisotropyRatio]);
    
    if INV.useSurfLapl == 1
        fprintf(fh,'%s\n','SL regularization');
    else
        fprintf(fh,'%s\n','inverse 1.distance^2 regularization ');
    end
    
    fprintf(fh,'mudep: %0.4g    murep: %0.4g \n',OPT.DEP.mu,OPT.REP.mu);
    fprintf(fh,'%s\n','ifix(1) tfix(1) dampd   dampr usetimes   tblc   mode  maxiter');
end

% Construct empty matrix to log results:
RESNOW = [];

%% PRELUDE TO ITERATIVE OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute reg and res and tres for initial estimate
if INV.doSave
    fprintf(fh,'%s\n','iter   min    mean     max     STD    reg     rd     tresd');
    temp = [0 min(OPT.DEP.tims) diff(range(OPT.DEP.tims)) max(OPT.DEP.tims) std(OPT.DEP.tims) TST.reg TST.rd TST.tres];
    sprintf(   '%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f',temp(1:8));
    fprintf(fh,'%3d %7.1f %7.1f %7.1f  %6.1f %6.1f %7.4f %7.4f\n',temp(1:8));
    fprintf(fh,'%s\n','iter   min   range   max    std  minAPD rangeAPD maxAPD stdAPD  reg   rd     tres');
end

disp('iter   min   range  max    STD  minAPD maxAPD stdAPD   reg      rd     tresd');

if isempty(diff(range((OPT.DEP.tims)))), dro = NaN; else; dro = diff(range((OPT.DEP.tims))); end

temp = [0 min(OPT.DEP.tims) dro max(OPT.DEP.tims) std(OPT.DEP.tims)...
    min(OPT.REP.tims-OPT.DEP.tims) ...
    max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims)...
    TST.reg TST.rd TST.tres];
fprintf('%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f Start time: %s \n',temp,datestr(now,'HH,MM.SS'));
% disp(txt);

if INV.doSave
    fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.3f %6.3f',temp);
    fprintf(fh,'time run: %s\n',datestr(now,'HH,MM.SS')) ;
end

RESNOW(1,:) = [temp norm(TST.RES)]; % File results of initial estimation in first row of log-matrix:
ik          = 2;                    % set count for iteration step to 2.

%% iterative approach to solving non-linear parameter estimation alternate between dep and rep
iter        = 0;                                                                                        % set count of iterations to zero
score       = 1;                                                                                        % start outer loop
startTime   = clock;                                                                                    % mark time of start inverse calculation
depscore    = 1;                                                                                        % marker to indicate if a next iteration step for depolarization is desired
repscore    = 1;                                                                                        % marker to indicate if a next iteration step for repolarization is desired

if INV.mode == 1, loopstart = 1; nloop = 1; repscore = 0; end                                           % dep only
if INV.mode == 2, loopstart = 2; nloop = 2; depscore = 0; end                                           % rep only
if INV.mode == 3, loopstart = 3; nloop = 3; end                                                         % notch only
if INV.mode == 4, loopstart = 1; nloop = 2; end                                                         % dep & rep
if INV.mode == 5, loopstart = 1; nloop = 2; end                                                         % dep & notch & rep
if INV.mode == 6, loopstart = 1; nloop = 1; repscore = 0; end                                           % dep & notch
if INV.mode == 7, loopstart = 7; nloop = 7; end                                                         % AP amplitude only

totscore    = 0;                                                                                        % Set total score of succesfull iteration steps to zero
% TST.rd      = 1;                                                                                        % Set initial RD value to 1
% Turned off - JvdW 9-5-2018

while iter < miniter || iter < INV.maxiter && (depscore || repscore) &&  TST.rd > INV.MINRD             % Check if a new iteration stepped is required
    iter = iter + 1;                                                                                    % increase count of iterations with one
    for loop = loopstart:nloop                                                                          % Initiate loop for mode's define above in INV.mode
        switch loop
            case 1
                [depscore,OPT.DEP,TST] = optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);           % Try to find depolarization times that better fit with the measured signal
                if INV.mode == 5 || INV.mode == 6
                    if INV.doAmplitude
                        [tmpscore,OPT.AMP,TST] = optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);     % Amplitude
                        warning('MATLAB:UndefinedFunction','optimizeAmpl.m is not fully checked!!')
                    else
                        [tmpscore,OPT.NOT,TST] = optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);    % Notch
                        warning('MATLAB:UndefinedFunction','optimizeNotch.m is not fully checked!!')
                    end
                end
                
                opt                     = OPT.DEP;                                                      % set [opt] to new depolarization results
                score                   = depscore;                                                     % File depscore to see if iterations should continue
                totscore                = totscore + depscore;                                        	% Increase total score with depscore
                
            case 2
                if INV.mode == 8
                    [tmpscore,OPT.NOT]  = optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);           % Notch
                    warning('MATLAB:UndefinedFunction','optimizeNotch.m is not fully checked!!')
                end
                [repscore,OPT.REP,TST]  = optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);          % Repolarization
                
                opt                     = OPT.REP;                                                      % set [opt] to new repolarization results
                score                   = repscore;                                                     % File repscore to see if iterations should continue
                totscore                = totscore + repscore;                                        	% Increase total score with repscore
                
            case 3
                [tmpscore,OPT.NOT,TST] = optimizeNotch(INV,OPT.NOT,OPT.DEP,OPT.REP,OPT.AMP);            % Notch
                warning('MATLAB:UndefinedFunction','optimizeNotch.m is not fully checked!!')
                opt = OPT.DEP;                                                                          % set [opt] to new depolarization results
                
            case 7
                [depscore,OPT.DEP]      = optimizeDepRep(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);          % Depolarization
                [tmpscore,OPT.AMP]      = optimizeAmpl(INV,OPT.AMP,OPT.DEP,OPT.REP,OPT.NOT);            % Amplitude
                warning('MATLAB:UndefinedFunction','optimizeAmpl.m is not fully checked!!')
                [repscore,OPT.REP,TST]  = optimizeDepRep(INV,OPT.REP,OPT.DEP,OPT.NOT,OPT.AMP);          % Repolarization
                
                opt                     = OPT.DEP;                                                      % set [opt] to new depolarization results
                score                   = depscore;                                                     % File depscore to see if iterations should continue
                totscore                = totscore + depscore;                                        	% Increase total score with depscore
        end
        
        if isempty(diff(range(opt.tims))), dro = NaN; else; dro = diff(range(opt.tims)); end             %
        
        temp = [iter min(opt.tims) dro max(opt.tims) std(opt.tims) min(OPT.REP.tims-OPT.DEP.tims) ...   % put all results in temp vector
            max(OPT.REP.tims-OPT.DEP.tims) std(OPT.REP.tims-OPT.DEP.tims) TST.reg TST.rd TST.tres];
        
        if score == 0                                                                                  % Check if results increased with last iteration step, if not change temp:
            k           = rem(ik,2)+1;
            temp(9:10)  = RESNOW(ik-k,9:10);
            temp(11)    = RESNOW(ik-1,11);
        end
        
        % Show results on screen:
        txt = sprintf('%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.2f %6.3f %6.3f run %s time %s',...
            temp,opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
        disp(txt);
        
        if INV.doSave                                                                                 % If doSave = 1, write results to ASCI-file
            fprintf(fh,'%3d %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.2f %6.3f %6.3f',temp);
            fprintf(fh,'run %s time %s\n',opt.type,datestr(datenum(clock)-datenum(startTime),'HH,MM.SS'));
        end
        
        RESNOW(ik,:)    = [temp norm(TST.RES)];                                                         % File results of last iteration step
        ik              = ik+1;                                                                         % Increase count for iteration with 1
    end
end

fprintf('lambda: %f\n',OPT.DEP.lambopt);

%%%% end outer loop; (loop1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Put final results into meas-structure:
TST             = gettres_v(INV,OPT.DEP,OPT.REP,OPT.NOT,OPT.AMP);                                       % Determine how good the fit with final depolarization and repolarization times is.
meas.depfinal   = OPT.DEP.tims;                                                                         % Final depolarization times
meas.repfinal   = OPT.REP.tims;                                                                         % Final repolarization times

if INV.useAmpl,     meas.amplfinal      = OPT.AMP.pot; end                                              % If Ampl is used: file amplitude potential values
if INV.useNotch,    meas.notchPotfinal  = OPT.NOT.pot; end                                              % If notch is used: file notch potential values

meas.rdfinal    = TST.rd;                                                                               % Final RD value
meas.regfinal   = TST.reg;                                                                              % Final regularization term
meas.tresfinal  = TST.tres;                                                                             % Final residue of RD plus regularization
COR             = corrcoef(TST.PHIA,INV.PHIREF);                                                        % Calculate correlation
meas.corfinal   = COR(2,1);                                                                             % Final correlation value between simulated and measured signals
meas.iterfinal  = iter;                                                                                 % Final number of iterations
meas.log        = RESNOW;                                                                               % Log-file

if INV.doSave               % If doSave = 1, write results to ASCI-file
    fprintf(fh,'\n Performance \n');
    fprintf(fh,'rdfinal:  %1.4f \n',meas.rdfinal);
    fprintf(fh,'corfinal: %1.4f \n',meas.corfinal);
    saveasci(outfil,[meas.depfinal meas.repfinal meas.repfinal-meas.depfinal]); 
    fclose(fh); 
end                          
%% =======================================================================