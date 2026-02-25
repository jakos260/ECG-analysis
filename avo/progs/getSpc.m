function S = getSpc(T,dep,rep,SPECS,mode,ADJneigh)

% S =getSp(dep,rep,size(T,2),ADJneigh,mode, SPECS)';
if mode == 1
    S =getSp(1,size(T,2),ADJneigh, dep,rep,[SPECS.initialSlope SPECS.plateauslope SPECS.repslope] )';
else
    S =getSp(3,size(T,2),ADJneigh, dep,rep,[SPECS.initialSlope SPECS.plateauslope SPECS.repslope] )';
end