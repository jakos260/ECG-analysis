function GEOM = selectLeads(GEOM,useLeads,useZeromean)

% GEOM.LAY = [ 3 3 0;...
%            1 1 1; 1 2 2; 1 3 3;...
%            2 1 4; 2 2 5; 2 3 6;...
%            3 1 7; 3 2 8;3 3 9];

GEOM.BSMall = GEOM.BSM;
GEOM.AMAall = GEOM.AMA;
if max(useLeads) > size(GEOM.BSM,1)
    useBSMleads = 1:size(GEOM.BSM,1);
else
    useBSMleads= useLeads;
end
if useZeromean==1
    GEOM.BSM = zeromean(GEOM.BSM(useBSMleads,:));
    GEOM.AMA = zeromean(GEOM.AMA(useLeads,:));
elseif useZeromean==-2 %assumed measured wct    
    GEOM.BSM    = GEOM.BSM(useBSMleads,:);
    GEOM.AMA    = GEOM.AMA(useLeads,:) - ones(length(useLeads),1)*mean(GEOM.AMA(1:3,:));
    GEOM.AMAall = GEOM.AMAall - ones(size(GEOM.AMAall,1),1)*mean(GEOM.AMA(1:3,:));
elseif useZeromean==2 %wct    
    GEOM.BSM    = GEOM.BSM(useBSMleads,:) - ones(length(useBSMleads),1)*mean(GEOM.BSM(1:3,:));
    GEOM.AMA    = GEOM.AMA(useLeads,:) - ones(length(useLeads),1)*mean(GEOM.AMA(1:3,:));
    GEOM.AMAall = GEOM.AMAall - ones(size(GEOM.AMAall,1),1)*mean(GEOM.AMA(1:3,:));
    % bipolar leads I II III
%     GEOM.BSM    = [GEOM.BSM; GEOM.BSM(2,:) - GEOM.BSM(1,:); GEOM.BSM(3,:) - GEOM.BSM(1,:); GEOM.BSM(3,:) - GEOM.BSM(2,:)];
%     GEOM.AMA    = [GEOM.AMA; GEOM.AMA(2,:) - GEOM.AMA(1,:); GEOM.AMA(3,:) - GEOM.AMA(1,:); GEOM.AMA(3,:) - GEOM.AMA(2,:)];
else
    GEOM.BSM = GEOM.BSM(useBSMleads,:);
    GEOM.AMA = GEOM.AMA(useLeads,:);
end

remove = 1:size(GEOM.BSMall,1);
remove(useLeads)=[];

if size(GEOM.LAY,2) > length(useLeads)
    LAY = [GEOM.LAY(1,:); GEOM.LAY(length(useLeads)+1,:)];
    for i=length(remove):-1:1
        LAY(LAY(:,3) >remove(i),3) = LAY(LAY(:,3) >remove(i),3) - 1;
    end

    GEOM.LAY=LAY;
end

if strcmp(GEOM.type,'ventricles')
else
end