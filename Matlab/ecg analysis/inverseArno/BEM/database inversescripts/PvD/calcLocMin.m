function locmin=calcLocMin(DIST,dep,locmin)

locMinArea=30;
if nargin<3 || isempty(locmin)
    locmin=[];
    for i=1:length(dep)
        if all(dep(DIST(:,i)<locMinArea)-dep(i)>=0)		
            locmin=[locmin i];
        end
    end
else
    
    locs=[];
    for i=1:length(locmin)
        locs = [locs find(DIST(locmin(i),:)< locMinArea)];
    end
    locs =unique(locs);
    orglocmin = locmin ;
    locmin=[];
    for i=1:length(locs)
        if all(dep(DIST(:,locs(i))<locMinArea) - dep(locs(i))>=0)		
            locmin=[locmin locs(i)];
        end
    end
    if length(orglocmin) ~= length(locmin)
        a=1;
    end
end