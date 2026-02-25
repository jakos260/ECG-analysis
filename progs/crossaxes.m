% crossaxes.m
% call to crossub should preceed;
% used to set desired axes and labels to cross-section plot
% 20062019

   mark_mean=plot(mean(VER(:,2)),mean(VER(:,1)),'r+');
   xlabel('y-coordinates')
   ylabel('x-coordinates')
   
   fixwin=get(ie7,'val');
   if fixwin==1,
      axis normal
   end
   
   labs=get(gca,'Yticklabel');
   [nlabs jlabs]=size(labs);
   newlabs=labs;
   for i=1:nlabs,
       newlab=sscanf(labs(i,:),'%f');
       lll=sprintf('%f',-newlab);
       newlabs(i,:)=lll(1:jlabs);
    end
    set(gca,'Yticklabel',newlabs)
    fixedassen=axis;

tit=title(['file: ' geom],'interpreter','none');
set(ie1,'string',num2str(zlevel))