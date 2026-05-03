function [avgshift shift]=alignsignals(X,Y,maxlag,plotflag)
if size(X,1)~=size(Y,1) % Assuming rows are channels, columns are samples
    if size(X,2)~=size(Y,2)
        error('Different number of channels');
    else
        X=X';
        Y=Y';
    end
end

if ~exist('plotflag','var')
    plotflag=0;
end
if ~exist('maxlag' )
    maxlag=[]; % empty maxlag is ignored by xcorr
end
if plotflag, figure; end
shift=zeros(1,size(X,1));

for i=1:size(X,1)
    xc=xcorr(X(i,:),Y(i,:),maxlag);
    t=(1:length(xc))-(length(xc)+1)/2;
    [xcmax imax]=max(xc);
    shift(i)=t(imax);
    if plotflag
        subplot(2,1,1);
        hold off
        plot(X(i,:),'r');
        hold on
        plot(Y(i,:),'b');
        legend('X','Y');
        subplot(2,1,2);
        hold off
        plot(t,xc);
        hold on
        plot(shift(i),xcmax,'*r');
        text(shift(i),xcmax,num2str(shift(i),'%0.1f'));
        pause (0.5);
    end
end
avgshift=mean(shift);

