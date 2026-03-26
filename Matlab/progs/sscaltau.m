% sscaltau.m
% shift and scale tau values to minimize rd over the interval tmin tmax
  %dep=timsd; rep=timsr;
  
  pol='all';
  RD=[];
  normph=norm(PHI(:,tmin:tmax),'fro');
  for ii=15:19;
      shift=ii-2;
      for jj=1:5;
          scal=1.1^(jj-7);
          dep=shift+scal*timsd;
          dep=max(dep,1);
          rep=ttm-0.4*(dep-mean(dep));
          gets;
          PHIA=A*S;
          RD=[RD; ii jj shift scal norm(PHI(:,tmin:tmax)-PHIA(:,tmin:tmax),'fro')/normph];
      end
  end
 RD
  