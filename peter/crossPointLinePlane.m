function [crosspoint, lambda, crossed]= crossPointLinePlane(planevector,planePosIn,v1,v2)

lambda = zeros(size(planePosIn,1),1);
crossed = zeros(size(planePosIn,1),1);
crosspoint = zeros(size(planePosIn,1),3);

toremove=[];
for i=1:size(planePosIn,1)
    planePos = planePosIn(i,:);
    if planePos == v1
        crosspoint(i,:) = v1;
        lambda(i) = 0;
        crossed(i) = 1;
    elseif planePos == v1
        crosspoint(i,:) = v2;
        lambda(i) = 1;
        crossed(i) = 1;
    else
        dl= v2 - v1;
        denom = dot(planevector,dl);
        numtop = dot(planevector,planePos-v1);
        if abs(denom) > eps && ... % NOT parallel no cross-section
           abs(numtop) > eps       % NOT parallel in plane
            crossed(i) = 1;
            lambda(i) = numtop / denom;
            if  lambda(i) < 0 || lambda(i) > 1
                crossed(i) = 2;
                toremove=[toremove;i];
            end
            crosspoint(i,:) = v1 + dl * lambda(i);           
        else
            toremove=[toremove;i];
        end
    end
end

crossed(toremove)=[];
crosspoint(toremove,:)=[];
lambda(toremove)=[];