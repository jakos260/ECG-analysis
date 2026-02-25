% tri_cleaner_plus.m
% 20140409
% function [VER,ITRI,list]=tri_cleaner_plus(VERIN,ITRIIN)
% PLUS: first search duplicate vertices; and remove duplicates

% remove any stray vertices from VER;
% PLUS: spot and remove duplicate vertices
% relabel ITRI accordingly

function [VER,ITRI]=tri_cleaner_plus(VER,ITRI)
% search duplicate VERTICES

nver=size(VER,1);
ntris=size(ITRI,1);
DUPLO=[];
tol=1.e-5;
for i=1:nver,
    for j=i+1:nver,
        if norm(VER(i,:)-VER(j,:))<tol,
            DUPLO=[DUPLO; [i j]];
        end
    end
end
nduplos=size(DUPLO,1);

if ~isempty(DUPLO),
    ivers=ITRI(:); % all vertex labels in ITRI
    for i=1:nduplos,
        k=find(ivers==DUPLO(i,2));
        if ~isempty(k),
            ivers(k)=DUPLO(i,1);
        end
    end
    ITRI=reshape(ivers,ntris,3);
end

  vers=ITRI(:);
  vers=unique(vers); % all vertices named in ITRI
  nverintri=size(vers,1);
  nverin=size(VER,1);
  if nverin==nverintri, return,end
  
  LIST=[(1:nverin)' zeros(nverin,1)];
  LIST(vers,2)=1;
  VERLIST=[LIST(LIST(:,2)==1,1) (1:nverintri)'];
  VER=VER(VERLIST(:,1),:);
  INDEXLIST=[(1:nverin)' zeros(nverin,1)];
  INDEXLIST(VERLIST(:,1),2)=INDEXLIST(VERLIST(:,2));

  ITRI=[INDEXLIST(ITRI(:,1),2) INDEXLIST(ITRI(:,2),2) INDEXLIST(ITRI(:,3),2)]; 
  %list=VERLIST(:,1);


