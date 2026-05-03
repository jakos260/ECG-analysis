
function meas=inverseQuasiNewton(varargin)%GEOM,inittimes,initcase,runcase,usetimes,leads,mode,mudep,murep)
% inverse determination of timing of Equivalent Double Layer source EDL;
% optimization based on basic Marquardt type of solving thenon-linear
% estimation problem alternate iterations between depolarization and
% repolarization
% Peter van Dam; 2010 november.
% All rights reserved Peacs, Arnhem  the Netherlands

% if INV.mode==1, loopstart=1;nloop=1; repscore=0;end % dep only
% if INV.mode==2, loopstart=2;nloop=2; depscore=0;end % rep only
% if INV.mode==3, loopstart=3;nloop=3; end % notch
% if INV.mode==4, loopstart=1;nloop=2; end % dep rep
% if INV.mode==5, loopstart=1;nloop=3; end % dep rep & notch
% if INV.mode==6, loopstart=6;nloop=6; end % AP amplitude only

% inverse procedure parameters
INV.mudep           = 0.01;
INV.murep           = 1.5e-4;
INV.muampl          = 1.5e-4;

INV.repOpt          = 'apd';
INV.MINRD           = 0.15;
INV.maxiter         = 25;
INV.lambopt         = 0.1;
INV.useAmpl         = 0;
INV.useNotch        = 0;
INV.useSurfLapl     = 1;
INV.useWeighedRd    = 1;
INV.experimentpvd   = 0;
INV.casedir         = 'invresults\';

INV.subname         = '';
INV.leads           = 1:65;
INV.mode            = 4;
INV.doWeight        = 0;
INV.doAmplitude     = 0;
INV.reg             = [];
INV.blmode=1; % oostep1, 1 blcorrection qrsonset-T, 2: (Vstim-10)-T.
INV.lpass=5;
INV.useAntiHakkel=false;
INV.geomdir='';

% criteria for ending iterations
INV.lambdamax=500;
INV.stopcrit=2e-4;

%%
if length(varargin) < 3
    error('This routine needs at least three parameters');
else
    GEOM = varargin{1};
    INV.amplitude = ones(size(GEOM.VER,1),1);
    INV.initdep = varargin{2};
    INV.initrep = varargin{3};
    INV.leads = 1:length(GEOM.TVER);
    usetimes = size(GEOM.BSM,2) - GEOM.SPECS.onsetqrs;
    pp=4;
    while pp<=nargin
        if ischar(varargin{pp})
            key=lower(varargin{pp});
            switch key
                case 'mudep'
                    INV.mudep = varargin{pp+1};pp=pp+2;
                case 'murep'
                    INV.murep = varargin{pp+1};pp=pp+2;
                case 'repopt'
                    INV.repOpt = varargin{pp+1};pp=pp+2;
                case 'muampl'
                    INV.muampl = varargin{pp+1};pp=pp+2;
                case 'estimateampl'
                    INV.useAmpl = varargin{pp+1};pp=pp+2;
                case 'estimatedamplitude'
                    amplitude = varargin{pp+1};pp=pp+2;
                case 'estimatenotch'
                    useNotch = varargin{pp+1};pp=pp+2;
                case 'maxiter'
                    INV.maxiter = varargin{pp+1};pp=pp+2;
                case 'minrd'
                    INV.MINRD = varargin{pp+1};pp=pp+2;
                case 'casedir'
                    INV.casedir = varargin{pp+1};pp=pp+2;
                case 'logname'
                    INV.subname = varargin{pp+1};pp=pp+2;
                case 'usetimes'
                    INV.usetimes = varargin{pp+1};pp=pp+2;
                case 'leads'
                    INV.leads = varargin{pp+1};pp=pp+2;
                    %                 case 'weighed'
                    %                     INV.doWeight = varargin{pp+1};pp=pp+2;
                case 'weighedrd'
                    INV.useWeighedRd  = varargin{pp+1};pp=pp+2;
                case 'mode'
                    INV.mode = varargin{pp+1};pp=pp+2;
                case 'reg'
                    INV.reg = varargin{pp+1};pp=pp+2;
                case 'blmode'
                    INV.blmode = varargin{pp+1};pp=pp+2;
                case 'lpass'
                    INV.lpass = varargin{pp+1};pp=pp+2;
                case 'antihakkel'
                    INV.useAntiHakkel = varargin{pp+1};pp=pp+2;
                case 'ctcpath'
                    INV.ctcpath = varargin{pp+1};pp=pp+2;
                    
                otherwise
                    warning('unknown parameter%s',key);
                    pp=pp+2;
            end
        end
    end
end

%% Inverse Quasi-Newton:

qnpath = '/Users/arnojanssen/Documents/stw/BEM/inverse/qn/bin';

delete(fullfile(qnpath,'invml*'));
copyfile(INV.ctcpath, fullfile(qnpath, 'invml.ctc'));
savematold(fullfile(qnpath, 'invml.ama'), GEOM.AMA);          % Amplitude??
saveasci(fullfile(qnpath,'invml.est'), INV.initdep);
savetri(fullfile(qnpath,'invml.tri'), GEOM.VER, GEOM.ITRI);

if strcmp(GEOM.type,'_atria') % no baseline correction for the atria
    
    if INV.blmode ~= 1, warning('Only baselinemode 1 supported for atria'); end
    
    INV.BSM =  baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetP:GEOM.SPECS.endtwave), INV.lpass));
    
    if INV.mode == 1
        INV.usetimes = min(GEOM.SPECS.onsetqrs,size(INV.BSM,2));
    else
        INV.usetimes = min(usetimes,size(INV.BSM,2));
    end
    
    INV.PHIREF = INV.BSM(:,1:INV.usetimes);

else
    
    usetimes = size(GEOM.BSM,2) - GEOM.SPECS.onsetqrs;
    
    if INV.blmode == 0,
        INV.BSM = lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), INV.lpass);
    elseif INV.blmode == 1,
        INV.BSM =  baselinecor(lowpassma(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), INV.lpass));
    elseif INV.blmode == 2,
        bsmt    = baselinecor(GEOM.BSM,GEOM.SPECS.time_Vstim-10,GEOM.SPECS.endtwave);
        INV.BSM = lowpassma(bsmt(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave),INV.lpass);
    end
    
    if INV.mode == 1,
        usetimes = min(GEOM.SPECS.qrsduration,size(INV.BSM,2));
    else
        usetimes = min(usetimes,size(INV.BSM,2));
    end
    
end

PHIREF = INV.BSM(:,1:usetimes);

savematold(fullfile(qnpath,'invmlecg.mes'), PHIREF);

oldpath = pwd;
cd(qnpath);

inppath = fullfile(qnpath,'invml.inp');
fp      = fopen(inppath,'w');

fprintf(fp,'%f\n%d\n%d\n%d\ny\ny\ninvml.est\n*\ninvml.ama\ninvmlecg.mes\ninvml.ctc\ninvml.tri\ninvml.mon\ninvml.tim\n', INV.mudep,INV.maxiter,INV.lpass,usetimes);

fclose(fp);
pause(1);

system(sprintf('%s < %s >>NULL',fullfile(qnpath,'qn') ,inppath));

meas.depfinal   = loadmat('invml.tim');
meas.PSIA       = loadmat('invml.tim.sim');
fp              = fopen('invml.mon','r');
meas.log        = fread(fp,'char');

fclose(fp);
cd(oldpath);

meas.repfinal = meas.depfinal+400;


return

