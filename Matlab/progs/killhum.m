function CLEANX = killhum(X,f,fsample,mu);

% CLEANX = killhum(X,f,fsample,mu);
% adaptive filter for harmonic component suppression; notch type filter
% X:      signal matrix (nsig,nt) 
% f:      notchfrequency; fsample: samplingfrequencey 
% mu:     adaptation coef. (< 1, typical 0.01 - 0.1);
%         notchwidth: 2*mu*f
%         small values of mu reduce the width of the notch at the expense of a
%         longer adaptation time; empirically: tau=pi/(2*mu) 
% CLEANX: the filtered output;
% vectorized version for nsig simultaneously recorded signals; 

% B. Widrow & S.D. Sterns: Adaptive signal Processing; Prentice Hall; 1985; 
% chapt12; page 317;
% algorithm sequentially fits a linear combination of cosine and sine at notch
% freqency to the input, and outputs the misfit = residual difference.

% A. van Oosterom; 21-9-2004

[nsig nt]=size(X);

CLEANX = zeros(nsig,nt); 
dfi=2*pi*f/fsample;
W=zeros(nsig,2);% assumes no hum to be present
for k=1:nt,
    y=[cos(k*dfi);sin(k*dfi)];
    hum = W*y;
    CLEANX(:,k) = X(:,k) - hum;
    DELW=mu*CLEANX(:,k)*y';
    W=W+DELW;
end

% scaling correction 
CLEANX=CLEANX*(1-mu/2);

