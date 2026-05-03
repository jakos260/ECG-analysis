% function [n,x]=histo(vals,fig,xbins,col,linewidth)
% draw histogram on figure(fig) ; normalized counts
% xbins should specify centers of the desired bins, having equal width
% if nargin<3, xbins=[min(vals):(max(vals)-min(vals))/11:max(vals)]; end
% A. van Oosterom
% 20121103

function [N,X]=histo(vals,fig,xbins,col,linewidth)
figure(fig) 
if nargin<5, linewidth=1; end
if nargin<4, col='b'; end
if nargin<3, xbins=[min(vals):(max(vals)-min(vals))/11:max(vals)]; end

[N, X] = hist(vals,xbins);
N=N/max(size(vals)*(xbins(2)-xbins(1)));
h = bar(X, N, 'hist');
set(h,'Facecolor','none','Edgecolor',col,'LineWidth',linewidth)
xlabel('x')
ylabel('fraction/binsize')
