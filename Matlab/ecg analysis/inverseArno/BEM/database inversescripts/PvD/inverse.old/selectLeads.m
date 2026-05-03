function GEOM = selectLeads(GEOM,useLeads,useZeromean)

% GEOM.LAY = [ 3 3 0;...
%            1 1 1; 1 2 2; 1 3 3;...
%            2 1 4; 2 2 5; 2 3 6;...
%            3 1 7; 3 2 8;3 3 9];

GEOM.BSMall = GEOM.BSM;
GEOM.AMAall = GEOM.AMA;

if useZeromean
    GEOM.BSM = zeromean(GEOM.BSM(useLeads,:));
    GEOM.AMA = zeromean(GEOM.AMA(useLeads,:));
else
    GEOM.BSM = GEOM.BSM(useLeads,:);
    GEOM.AMA = GEOM.AMA(useLeads,:);
end