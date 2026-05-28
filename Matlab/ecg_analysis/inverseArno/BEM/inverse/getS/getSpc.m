% function S = getSpc(T,depIn,repIn,SPECS,mode,ADJneigh)
%
% This function simulates the Transmembrane potentials [TMP's] for the
% given depolarization & repolarization times, given the specified
% adjacency matrix and TMP spefications in SPECS.
%
% INPUT:
% T         = matrix of # vertices * timescale QRT
% depIn     = depolarization times [Has to be same size as ADJneigh, or length == 1]
% repIn     = repolarization times [Has to be same size as ADJneigh, or length == 1]
% SPECS     = Specs from prepareECG.m: Necessary are --> a) initial repslope, b) plateau slope and c) rep slope
% mode      = use only depolarization [1] or also repolarization [3]
% ADJneigh  = adjacency matrix with graphdist(ITRI)
%
% OUTPUT:
% S         = Simulated TMP's.
%
% Version 1: 01-apr-2015

function S = getSpc(T,depIn,repIn,SPECS,mode,ADJneigh)

% Check if length == 1:
if length(depIn) == 1, depIn = depIn * ones(length(ADJneigh),1); end
if length(repIn) == 1, repIn = repIn * ones(length(ADJneigh),1); end

% determine minimum depolarization time:
dep0 = min(depIn);

% check for NaN values:
if isnan(dep0), error('dep0 is NaN'); end

%% define input for getSp.c: 
if dep0 < 0,    extr_dep = dep0;    else; extr_dep = 0; end                                              % determine 

% mode 1: only depolarization, mode 3: S using cumsum (avo), requires a)
% initial repslope, b) plateau slope and c) rep slope
if mode < 0
    inp_mode = 2;               repparams   = [abs(SPECS.plateauslope) abs(SPECS.repslope)];            % parameters to describe TMP.
elseif mode == 1
    inp_mode = 1;               repparams   = [];                                                       % no parameters to describe TMP needed.
elseif SPECS.useCumsum == 1
    inp_mode = 3;               repparams   = [SPECS.initialSlope SPECS.plateauslope SPECS.repslope];   % parameters to describe TMP.
elseif SPECS.useCumsum == 0
    inp_mode = 4;               repparams   = [SPECS.initialSlope SPECS.plateauslope SPECS.repslope];   % parameters to describe TMP.
end

maxt        = size(T,2) - extr_dep;                                                                     % maximum time, is the number of columns of S. Corrected with earliest dep time.
dep_Sp      = depIn - extr_dep;                                                                         % Correct dep times with earliest dep time.
rep_Sp      = repIn - extr_dep;                                                                         % Correct rep times with earliest dep time.   

%% Call c++ program getSp.c to simulate the TMP's based on the depolarization and repolarization times:
S = getSp(inp_mode, maxt, ADJneigh, dep_Sp, rep_Sp, repparams)';

% correct S for minimum depolarization smaller than zero.
if dep0 < 0, S = S(:,end-size(T,2)+1:end); end

% Check if length == 1, if yes bring output also back to length 1;
if length(depIn) == 1, S = S(1,:); end