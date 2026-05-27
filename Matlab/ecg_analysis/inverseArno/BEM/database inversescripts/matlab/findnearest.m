function [pnearest, distp, trinearest, vernearest, distver,la,mu] = findnearest(point,VER, ITRI, surftype, surfselect, n)
% function [pnearest, distp, trinearest, vernearest, distver,la,mu] = findnearest(point,VER, ITRI, surftype, surfselect, n)
% find nearest n-points (only n==1 working) on ventricular geometry (VER, ITRI, surftype) of type surfselect (string). n<0: find
% point on surface~=surfselect. 
lamu=[];
surfnames={'Epi','LVEndo','RVEndo','LVEndo','RVEndo','RVEndo','LVEndo'}; % Simplified

if ~exist('n','var'), n = 1; end

if n < 0, surfselect = unique(surfnames(~strcmp(surfselect,surfnames))); n = abs(n); end

if n > 1, error('only n == 1 implemented'); end

% Surftype 1-7: Epi, LVEndo,RVEndo,MitralValve, TricuspidValve,RVOT,Aorta
if ~exist('surfselect','var'), 
    surfignore = []; 
else
    % surfnames={ 'Epi', 'LVEndo','RVEndo','MitralValve',
    % 'TricuspidValve','RVOT','Aorta'};
    isurf = [];
    if ~iscell(surfselect)
        surfselect = {surfselect};
    end
    
    for i = 1:length(surfselect)
        isurf = [isurf find(strcmp(surfnames,surfselect{i}))];
    end
    
    if isempty(isurf)
        error('Unknown surface string');
    else
        surfignore =~ ismember(surftype,isurf);
    end
    
end

[LA,MU,dist]    = ptridist(point,VER,ITRI);
triignore       = all(ismember(ITRI,find(surfignore)),2);
dist(triignore) = inf;

[distp,imind]   = min(abs(dist));
pnearest        = (1-LA(imind)-MU(imind))*VER(ITRI(imind,1),:)+ LA(imind)*VER(ITRI(imind,2),:)+MU(imind)*VER(ITRI(imind,3),:);
trinearest      = imind;

for j = 1:size(VER,1);
    ndist(j) = norm(VER(j,:)-point);
end

ndist(surfignore)       = inf;
[distver,vernearest]    = min(ndist);
la                      = LA(imind);
mu                      = MU(imind); % return only mu of selected point
