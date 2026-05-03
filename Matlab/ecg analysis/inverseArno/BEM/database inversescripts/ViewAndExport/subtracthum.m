function sigout=subtracthum(sigin)
% 20130626 oostep1: original version

wwidth=20;
wstart=1;
period=sigin(:,wstart:(wstart+wwidth));
% period=baselinecor(period); % Use only when signal contains large baseline drift, else it wil just amplify quantisation noise  
period=period(:,wstart:wstart+wwidth-1);
period=period-repmat(mean(period,2),1,size(period,2));
nrep=ceil(size(sigin,2)/wwidth); % TBD align with start phase if wstart~=1
wave=repmat(period,1,nrep);
wave=wave(:,1:size(sigin,2));
% wave=lowpassma(wave(:,1:size(sigin,2)),3);
sigout=sigin-wave;