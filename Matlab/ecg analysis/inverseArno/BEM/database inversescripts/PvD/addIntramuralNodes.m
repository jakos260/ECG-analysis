function [ADJ3D,DIST3D] = addIntramuralNodes(VER,ITRI,ADJ,endoVER)

% this script creates locally transmural nodes and adds the connections to
% the adjacency matrix.I have limited the number of connections as much as
% possible. 


% load the geometry data
ADJSURF=graphdist(ITRI,VER,4);
ADJ3D=ADJ;
ADJ3D(ADJSURF>0)=0;
ADJ3D(ADJ3D<7)=0;

% create an ischemic zone in which the propagation velocity is changed.
% Only add the connectiosn within the ischemicZoneRadius distance from the
% centerIschemia 
centerIschemia = 1194; % node number
ischemicZoneRadius = 30; % mm

cDist = DIST(centerIschemia,:);
ADJ3D( ADJ3D > 40 ) = 0;
ADJ3D(cDist > ischemicZoneRadius ,cDist> ischemicZoneRadius )=0;
for i=1:length(VER)
    if any(ADJ3D(i,:)>0)
        ADJ3D(i,ADJ3D(i,:) > 2.0 * min(ADJ3D(i,ADJ3D(i,:)>0)) ) = 0;
        ADJ3D(:,i) = ADJ3D(i,:);
    end
end
% number of splits in a transmural connection
nSplit = 2;
% number of to be added vertices
Nnewver = nSplit * length(nonzeros(ADJ3D))/2;

% reserve memory
VERNEW = [VER; zeros(Nnewver,3)];
ADJ3DNEW = zeros(length(VERNEW));
ADJ3DNEW(1:length(VER),1:length(VER)) = ADJ;
ADJ3DNEW( ADJ3DNEW > 40 ) = 0;

% create the new vertices and its connections on any of th eselected
% transmural connections
verK=length(VER)+1;
for i=1:length(VER)
    for j=i:length(VER)
        if ADJ3D(i,j) > 0
            ADJ3DNEW(i,j) = 0;
            ADJ3DNEW(j,i) = 0;            
            % create the new vertices
            for k= 1:nSplit 
               VERNEW(verK,:) = VER(i,:) + (VER(j,:) - VER(i,:)) * k / (nSplit + 1);
               if k == 1
                   ADJ3DNEW(i,verK) = norm3d(VER(i,:) - VERNEW(verK,:));
                   ADJ3DNEW(verK,i) = ADJ3DNEW(i,verK);            
               elseif k==nSplit
                   ADJ3DNEW(j,verK) = norm3d(VER(j,:) - VERNEW(verK,:));
                   ADJ3DNEW(verK,j) = ADJ3DNEW(j,verK);                               
               else
                   ADJ3DNEW(verK-1,verK) = norm3d(VERNEW(verK,:) - VERNEW(verK-1,:));
                   ADJ3DNEW(verK,verK-1) = ADJ3DNEW(verK-1,verK);                                                 
               end
               verK = verK + 1;
            end
        end
    end
end

% create connections between the intramural connections (only for vertices 
% in close range to each other. This will introduce some numerical error,
% but probably limited for the purpose.
for i=1:length(VERNEW)
    for j=i:length(VERNEW)
        if ADJ3DNEW(i,j) == 0 
            d = norm3d(VERNEW(i,:) - VERNEW(j,:)); %mm
            if d < 5
                ADJ3DNEW(i,j) = d;
                ADJ3DNEW(j,i) = d;
            end
        end
    end
end


% create distance matrix for the intramural nodes too. time consuming!!!
DIST = graphdist(ADJ3DNEW); 