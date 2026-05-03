% function plotecg(GEOM, data, figname, scanmode, lpass, sigplotscale)
%
% Two figures are created and saved;
%   1. Final QRST: sigplot.m figure with measured [blue] and final inverse [red] QRST.
%   2. QRS complex: sigplot.m figure with measured [blue], initial
%      estimation [black] and final inverse solution [red] of QRS complex.
%
% INPUT:
% GEOM          = Structure created with InvInit.m. ; contains .SPECS, AMA & LAY 
% data          = Contains T-vector, initial estimation and final result of inverse calculation
% figname       = Name for figures to be saved as, inlcuding pathname. Without .jpg
% scanmode      = mode for getSmode [default == 1]
% lpass         = lowpass filtering parameter [lowpassma.m] [default == 1]
% sigplotscale  = parameter fro sigplot.m
%
% Version 1: 01-apr-2015

function plotecg(GEOM, data, figname, scanmode, lpass, sigplotscale)

if ~exist('scanmode','var'),        scanmode        = 1;    end
if ~exist('lpass','var'),           lpass           = 1;    end 
if ~exist('sigplotscale','var'),    sigplotscale 	= 1.5;  end 

sigplotleadnos  = 1; 	% Used for sigplot.m at the end of this script.

T               = data.T;
meas            = data.meas;
measinit        = data.measinit;

% Find init opposite:
PSIA    	= lowpassma(GEOM.AMA*getSmode(T, meas.depfinal, meas.repfinal, GEOM.SPECS, scanmode,GEOM.neigh), lpass);  % Construct final result simulated ECG
PSIAinit 	= lowpassma(GEOM.AMA*getSmode(T, measinit.dep, measinit.rep, GEOM.SPECS, scanmode, GEOM.neigh), lpass);   % Construct resulting simulated ECG from initial estimation

% Visualize final solution inverse QRST:
figure('Name','Final QRST');
clf;
sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.endtwave), ' ' , GEOM.LAY, sigplotscale, 'b', sigplotleadnos);
hold on
sigplot(PSIA,'blue = meas, red = sim', GEOM.LAY, sigplotscale, 'r', sigplotleadnos);

% save figure:
saveas(gcf, [figname '_QRST.jpg'])

% Visualize final solution inverse QRS:
figure('Name','Final QRS');
clf;
sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:(GEOM.SPECS.onsetqrs + GEOM.SPECS.time_Jpoint - 1)), '', GEOM.LAY, sigplotscale, 'b', sigplotleadnos);
hold on
sigplot(PSIA(:,1:GEOM.SPECS.time_Jpoint), '', GEOM.LAY, sigplotscale, 'r', sigplotleadnos);
hold on
sigplot(PSIAinit(:,1:GEOM.SPECS.time_Jpoint), 'blue = meas, red = sim, black = init', GEOM.LAY, sigplotscale, 'k', sigplotleadnos);

% Save figure:
saveas(gcf, [figname '_QRS.jpg'])