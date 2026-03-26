% time_run
secs=toc
if secs>=3600,
   hrs=floor(secs/3600);
   mns=floor(rem(secs,3600)/60);
   sec=secs-hrs*3600-mns*60;
   ['time run: '  num2str(hrs) ' hr; ' num2str(mns) ' min; ' num2str(sec) ' s']
else
   if secs>=60,
       mns=floor(secs/60);
       sec=secs-mns*60;
       ['time run: ' num2str(mns) ' min; ' num2str(sec) ' s']
   else
       ['time run: ' num2str(secs) ' s']
   end
end