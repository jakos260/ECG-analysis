

% function meas = qnopt(varargin)
% 
% GEOM
%
% 15-dec-2015

function meas = qnopt(varargin)

% inverse procedure parameters
INV.mudep           = 0.01;
INV.murep           = 1.5e-4;
INV.lpass           = 5;


% criteria for ending iterations
INV.lambdamax       = 500;
INV.stopcrit        = 2e-4;

%% Inverse Quasi-Newton:

qnpath = '/Users/arnojanssen/Documents/stw/BEM/inverse/qn/bin';

delete(fullfile(qnpath,'invml*'));                              % clear all old quasi-newton files
copyfile(INV.ctcpath, fullfile(qnpath, 'invml.ctc'));           % First create laplacian with lap-executable and then copy to this folder.
savematold(fullfile(qnpath, 'invml.ama'), GEOM.AMA);            % Save A-matrix
saveasci(fullfile(qnpath,'invml.est'), INV.initdep);            % Save initial estimation of depolarization values
savetri(fullfile(qnpath,'invml.tri'), GEOM.VER, GEOM.ITRI);     % Save ventricle surface mesh 
    

PHIREF = INV.BSM(:,1:usetimes);

savematold(fullfile(qnpath,'invmlecg.mes'), PHIREF);            % Save measured signal (to be used)

oldpath = pwd;
cd(qnpath);

inppath = fullfile(qnpath,'invml.inp');
fp      = fopen(inppath,'w');

fprintf(fp,'%f\n%d\n%d\n%d\ny\ny\ninvml.est\n*\ninvml.ama\ninvmlecg.mes\ninvml.ctc\ninvml.tri\ninvml.mon\ninvml.tim\n', INV.mudep, INV.maxiter, INV.lpass, usetimes);

fclose(fp);
pause(1);

system(sprintf('%s < %s >>NULL', fullfile(qnpath,'qn'), inppath));

meas.depfinal   = loadmat('invml.tim');                         % load final depolarization times
meas.PSIA       = loadmat('invml.tim.sim');                     % load corresponding reconstructed ECG [BSM]
fp              = fopen('invml.mon', 'r');                      % log-file of quasi-newton optimization
meas.log        = fread(fp, 'char');

fclose(fp);
cd(oldpath);

meas.repfinal = meas.depfinal + 400;
