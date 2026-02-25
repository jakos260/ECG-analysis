% function LOW=lowpassma(M,nwin)
% lowpass Moving Average filter over nwin>=2 samples, applied to rows of matrix M; finite impulse
% response; no phase shift
% zero output for frequencies (Hz) at multiples of  f_0 = 1/(nwin*sint)
% hence, nwin=1/(f_0*sint); sr=sampling_rate=1/sint.....  nwin=sr/f_0

% for highpass: use HIGH=M-lowpassma(M,nwin)

function LOW=lowpassma(M,nwin)
% 2005-02-17;  % V. Jacquemet/ A. van Oosterom
% 2008-02-12; AvO winb=round(... )   changed into:   winb=floor(....)

LOW=M;
if nwin <2, return, end
trans=0;
[nsig, nt]=size(M);
if nt==1,
    trans=1;
    M=M';
    [nsig, nt]=size(M);
end
winb=floor(nwin/2);
wine=nwin-winb;
LEAD=M(:,1)*ones(1,winb);
TRAIL=M(:,nt)*ones(1,wine);
M=[LEAD M TRAIL];
X = cumsum([zeros(nsig,1) M],2);
LOW = X(:,nwin+1:nwin+nt)-X(:,1:nt);
LOW = LOW / nwin;
if trans==1,
    LOW=LOW';
end