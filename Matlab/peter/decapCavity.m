function decapCavity(filedirIn)
    
% filedir = 'age65_w81_h170_new';
[VER,ITRI]=loadtri(fullfile(filedirIn,[filedirIn '_ventricles.tri']));
type=loadmat(fullfile(filedirIn,[filedirIn '_ventricles.typ']));
[LVER,LITRI]=loadtri(fullfile(filedirIn,[filedirIn '_lcav.tri']));
[RVER,RITRI]=loadtri(fullfile(filedirIn,[filedirIn '_rcav.tri']));
[TVER,TITRI]=loadtri(fullfile(filedirIn,[filedirIn '_thorax.tri']));

for i=length(LVER):-1:1
    if ~any(VER(:,1) == LVER(i,1) & ...
            VER(:,2) == LVER(i,2) & ...
            VER(:,3) == LVER(i,3) )
        [LVER,LITRI] = delNode(LVER,LITRI,i);
    end
end
[lver,litri]=getCap(VER,ITRI,type,[4 7]);litriOrg = litri;
for i=1:length(lver)
    ai = find(LVER(:,1) == lver(i,1) & ...
              LVER(:,2) == lver(i,2) & ...
              LVER(:,3) == lver(i,3) );
    if isempty(ai)
        LVER =[LVER; lver(i,:)];
        litri(litriOrg==i) = length(LVER);
    else
        litri(litriOrg==i) = ai;
    end
end

LITRI=[LITRI;litri(:,[2 1 3])];

for i=length(RVER):-1:1
    if ~any(VER(:,1) == RVER(i,1) & ...
            VER(:,2) == RVER(i,2) & ...
            VER(:,3) == RVER(i,3) )
        [RVER,RITRI] = delNode(RVER,RITRI,i);
    end
end
[rver,ritri]=getCap(VER,ITRI,type,5);ritriOrg=ritri;

for i=1:length(rver)
    ai = find(RVER(:,1) == rver(i,1) & ...
              RVER(:,2) == rver(i,2) & ...
              RVER(:,3) == rver(i,3) );
    if isempty(ai)
        RVER =[RVER; rver(i,:)];
        ritri(ritriOrg==i) = length(RVER);
    else
        ritri(ritriOrg==i) = ai;
    end
end
RITRI=[RITRI;ritri(:,[2 1 3])];

[pulmver,pulmitri]=getCap(VER,ITRI,type,6); pulmitriOrg = pulmitri;
for i=1:length(pulmver)
    ai = find(RVER(:,1) == pulmver(i,1) & ...
              RVER(:,2) == pulmver(i,2) & ...
              RVER(:,3) == pulmver(i,3) );
    if isempty(ai)
        RVER =[RVER; pulmver(i,:)];
        pulmitri(pulmitriOrg==i) = length(RVER);
    else
        pulmitri(pulmitriOrg==i) = ai;
    end
end
RITRI=[RITRI;pulmitri(:,[2 1 3])];

filedir = [filedirIn 'decap'];
mkdir(filedir);
savetri(fullfile(filedir,[filedir '_ventricles.tri']),VER,ITRI);
savemat(fullfile(filedir,[filedir '_ventricles.typ']),type);
savetri(fullfile(filedir,[filedir '_lcav.tri']),LVER,LITRI);
savetri(fullfile(filedir,[filedir '_rcav.tri']),RVER,RITRI);
savetri(fullfile(filedir,[filedir '_thorax.tri']),TVER,TITRI);




%%

function [capver,capitri]=getCap(VER,ITRI,type,ringType)
epiType=1;
if length(ringType)==2
    veri= find(type==ringType(1) | type==ringType(2) );
else
    veri= find(type==ringType);
end
adj=zeros(length(VER));
for i=1:length(ITRI)
    if length(ringType)==2
        a = find(type(ITRI(i,:)) == ringType(1) | type(ITRI(i,:)) == ringType(2)); 
    else
        a = find(type(ITRI(i,:)) == ringType);  
    end
    if length(a) == 2 && any(type(ITRI(i,:)) == epiType)
        if a(1) ==1 && a(2) == 3
            adj(ITRI(i,a(2)),ITRI(i,a(1)))=1;
        else
            adj(ITRI(i,a(1)),ITRI(i,a(2)))=1;
        end
    end
end
capver =VER(veri(1),:);
next=veri(1);
for i=2:length(veri)
    next = find(adj(next,:)==1);next=next(end);
    capver =[capver; VER(next,:)];    
end
capver = [capver; mean(capver)];
centeri = length(capver);
capitri=[centeri length(capver)-1 1 ];

for i=2:length(capver)-1
    capitri=[capitri;[centeri (i-1) (i) ]];
end
