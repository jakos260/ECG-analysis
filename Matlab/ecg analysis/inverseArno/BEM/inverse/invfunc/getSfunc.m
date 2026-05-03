% function S = getSfunc(T,dep,rep,SPECS,notchpot,scaleAmpl,mode,INV)
%
% Wrapper function for getSmode.m
%
% INPUT
% T         = matrix of (nr of vertices on ventricles mesh) x (duration BSM)
% dep       = depolarization times
% rep       = repolarization times
% SPECS     = INV.SPECS
% notchopt  = zeros, unless INV.useNotch = 1 in inversefunc.m
% scaleAmpl = Amplitudes, defined in inversefunc.m as [OPT.AMP.pot = INV.amplitude]
% mode      = mode of inverse calculation [1] dep, [4] rep etc.
% INV       = INV-structure constructed in inversefunc.m
%
% OUTPUT:
% S         = Simulated TMP's.
%
% Version 1: 01-apr-2015

function S = getSfunc(T,dep,rep,SPECS,notchpot,scaleAmpl,mode,INV)

% parameter:
notchparam = 50;

% Check if mode is set for repolarization [4] or any other mode:
if mode == 4,
    
    SPECS.scaleAmpl = scaleAmpl;                                    % Amplitudes
    
    if size(T,2) < max(rep) + notchparam,                           % Check if duration of ECG is smaller than maximum repolarization time + notchparam.
        Tt  = ones(length(dep),1)*(0:max(rep) + notchparam);        % If yes, adjust T;
        S   = getSmode(Tt,dep,rep,SPECS,mode);                    	% Get simulated TMP's in mode 4.
        S   = S(:,1:size(T,2));                                     % Adjust size of TMP's.
    else
        S   = getSmode(T,dep,rep,SPECS,mode,INV.neigh);             % Get simulated TMP's in mode 4.
    end
    
else
    SPECS.scaleAmpl = notchpot;                                     % Potentials at J-point?
    S               = getSmode(T,dep,rep,SPECS,mode,INV.neigh);     % Get simulated TMP's.
end