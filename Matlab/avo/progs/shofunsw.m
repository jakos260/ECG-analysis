% shofunsw.m
%call callback routine of triplot
funsw=get(ui1r4,'Value');
if funsw ==1,
  if exist('VALS')&isempty(VALS)==0,
     set(ui15,'vis','off')
     fun=VALS(:,column);
     set(hs,'FaceVertexCData',fun);
     set([ui16;ui17;ui3a;ui4],'vis','on'),
  else,
     set(ui15,'vis','on')
     set(ui1r4,'val',0);
     funsw=0;
  end
else, 
  fun=mean(VALS(:,column))*ones(nver,1);
  set(hs,'FaceVertexCData',fun);
  set([ui3a;ui4],'vis','off')
  set([ui15;ui16;ui17;],'vis','off')
end