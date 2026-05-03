% file trithinning.m
% date:0504018
      %G=tri2gra(VER,ITRI,1);
      %D=gra2dist(G);
      [A,D]=graphdist(VER,ITRI,0);
      maxd=max(max(D))
      
      node=1;
      iver=1:nver;
      select=node;
      % select second order neighbours of node
      for id=2:2:maxd,
         select=[select iver(D(node,:)==id)];
      end
      keep=select(1:2);
      select(1:2)=[];
      % remove any direct neighbours from select
   while isempty(select)==0,
       keep=[keep select(1)]
       bver=find(D(:,select(1))==1);
       select(ismember(select,bver))=[];
       select(1)=[];
   end
   text(VER(keep,1),VER(keep,2),VER(keep,3),'*','col','w');
% outputnodes in: keep