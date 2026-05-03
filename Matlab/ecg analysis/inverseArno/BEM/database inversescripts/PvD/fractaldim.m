function [df,sig]=fractaldim(varargin) 
%df=fractaldim(x,frame=60,resolutie=0.09,mode=1)

if length(varargin) < 1
	disp('No input signal!!!');
	return;
else
	frame=60;
	resolutie=0.09;
	mode=2;
	x									=varargin{1};
	if length(varargin) > 1, frame		=varargin{2}; end;
	if length(varargin) > 2, resolutie	=varargin{3}; end;
	if length(varargin) > 3, mode		=varargin{4}; end;
	if length(varargin) > 4, disp('no extra parameters required'); end
end


%=================== proposed method
L=1;
if mode==1
	df=zeros(size(x));
	mx1=ones(size(x));
	sig=diffrows(x');
	xx=sig./resolutie;
	xx=xx-min(xx);
	if round((max(xx)-min(xx)))+2 > 800
		error('too many boxes')
	end
	nul=zeros(round((max(xx)-min(xx)))+2);
	for i=frame+L:length(xx)
		dff=nul; 
		dff(round(xx(i-L-frame+1:i-L))+ 1 +...
			 round(xx(i-frame+1:i))*length(nul))=1;
		df(i)=sum(sum(dff));
	end
	df=100.*df./frame;
elseif mode==2
	df=zeros(size(x));
	mx1=zeros(size(x));	mx1(2:end)=x(1:end-1);
	sig=(x-mx1)./(x+mx1);
	xx=(ones(size(x))+sig)./resolutie;	
	if round((max(xx)))+2 > 800
		error('too many boxes')
	end
	nul=zeros(round((max(xx)))+2);
	for i=frame+L+2:length(x)
		dff=nul; 
		dff(round(xx(i-L-frame+1:i-L)) + 1 +...
			 round(xx(i-frame+1:i))*length(nul))=1;
		df(i)=sum(sum(dff));
	end	
	df=100.*df./frame;
else
	df=zeros(size(x));
	sig=x;
	xx=(x-min(x))./resolutie;
	if round((max(xx)-min(xx)))+1 > 800
		error('too many boxes')
	end

	nul=zeros(round((max(xx)-min(xx)))+1);
	for i=frame+L:length(x)
		dff=nul; 
		dff(round(xx(i-L-frame+1:i-L)) + 1 +...
			 round(xx(i-frame+1:i))*length(nul))=1;
		df(i)=sum(sum(dff));
	end	
	df=100.*df./frame;
end