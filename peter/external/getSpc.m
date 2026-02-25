function S = getSpc(T,depIn,repIn,SPECS,mode,ADJneigh)



% dep0= min(dep);
% if dep0 < 0
%     if mode < 0
%         S =getSp(2, size(T,2) - dep0,ADJneigh, dep-dep0, rep-dep0, [ abs(SPECS.plateauslope) abs(SPECS.repslope)] )';
%     elseif mode == 1
%         S =getSp(1,size(T,2) - dep0,ADJneigh, dep - dep0, rep-dep0, [SPECS.initialSlope SPECS.plateauslope SPECS.repslope])';
%     else
%         S =getSp(3, size(T,2) - dep0,ADJneigh, dep-dep0, rep-dep0, [SPECS.initialSlope SPECS.plateauslope SPECS.repslope] )';        
%     end
%     S= S(:,end-size(T,2)+1:end);
% else
    if length(depIn) == 1 
        dep = depIn * ones(length(ADJneigh),1);
    else
        dep=depIn;
    end
    if length(repIn) == 1 
        rep = repIn * ones(length(ADJneigh),1);
    else
        rep=repIn;
    end

    if mode < 0
        S =getSp(2, size(T,2),ADJneigh, dep, rep, [ abs(SPECS.plateauslope) abs(SPECS.repslope)] )';
    elseif mode == 1
        S =getSp(1,size(T,2),ADJneigh, dep )';
    elseif SPECS.useCumsum == 1
        S =getSp(3,size(T,2),ADJneigh, dep,rep,[SPECS.initialSlope SPECS.plateauslope SPECS.repslope] )';
    elseif SPECS.useCumsum == 0
        S =getSp(4,size(T,2),ADJneigh, dep,rep,[SPECS.initialSlope SPECS.plateauslope SPECS.repslope] )';
    end
    if length(depIn)==1
        S =S(1,:);
    end
% end