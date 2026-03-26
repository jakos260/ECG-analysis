% node2col.m
% script of triplot.m and newinit
% accept node as new column
column=node;

[idum maxcoll]=size(VALS);

if column<=maxcoll, 
   set(sltri,'value',column);
   displcol
end

next=0;
if next==1,
    
   set(sltri,'max',maxcoll');
   fun=VALS(:,column);
   set(hs,'FaceVertexCData',fun);
   hcbar=colorbar;
   hcbarpos=[.9 .12 .03 .6];
   set(hcbar,'position',hcbarpos)

   if symcolor==1,
      % symmetric color bar around zero
      extrmap=max(abs([fun;eps]));
      %startpos=get(gca,'Cameraposition');
      set(ax1a,'Clim',[-extrmap extrmap]);
      %hcbar=colorbar;
      set(hcbar,'position',hcbarpos),
   end
   
   if exist('ECGG'),
      % for initudl:
      figure(2);
      clf
      sigplot(ECG,' ',LAY,.7,'b',0);
      dep=TIMS(:,node);
      gets;
      PHI=AA*S(:,1:tmaxi);
      rd=norm(ECG-PHI,'fro')/norm(PHI,'fro');
      sigplot(PHI,' ',LAY,.7,'r',0), 

      % display rd
      ui35=uicontrol('style','text');
      uibox35=[.7 .9 .12 .05];
      set(ui35,'units','norm','position',uibox35,'vis','on','string',...
      [' rd=' num2str(rd)]);
   end 
   set(ui17,'string',num2str(node))


if exist('ui4'),
  set(ui4,'string',sprintf('%3.3f',fun(node)));
end
end