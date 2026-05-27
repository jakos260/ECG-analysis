function S = getSpc(T,dep,rep,SPECS,mode,ADJneigh)

dep0 = min(dep);
if isnan(dep0), error('dep0 is NaN'); end

if dep0 < 0
    if mode == 1
        S = getSp(1, size(T,2)-dep0, ADJneigh, dep-dep0, rep-dep0, [SPECS.initialSlope SPECS.plateauslope SPECS.repslope])';
    else
        S = getSp(3, size(T,2)-dep0, ADJneigh, dep-dep0, rep-dep0, [SPECS.initialSlope SPECS.plateauslope SPECS.repslope])';
    end
    S = S(:,end-size(T,2)+1:end);
else
    if mode == 1,
        S = getSp(1, size(T,2), ADJneigh, dep, rep, [SPECS.initialSlope SPECS.plateauslope SPECS.repslope])';
    else
        S = getSp(3, size(T,2), ADJneigh, dep, rep, [SPECS.initialSlope SPECS.plateauslope SPECS.repslope])';
    end
end