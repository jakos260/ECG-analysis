function [VER,ITRI]=delEdge(VER,ITRI,ver1,ver2)

if length(ver1)==1
   node1= ver1; 
else
    node1 = find(VER(:,1)==ver1(1) & VER(:,2)==ver1(2) & VER(:,3)==ver1(3) );
end
if length(ver2)==1
    nod2=ver2;
else
    node2 = find(VER(:,1)==ver2(1) & VER(:,2)==ver2(2) & VER(:,3)==ver2(3) );
end
newVer = ver1 + (ver2-ver1)/2;
VER(node1,:) = ones(length(node1),1)*newVer;
VER(node2,:) = ones(length(node2),1)*newVer;
return
if length(node1)==2 && length(node2)==2
    
    for i=1:length(node1)
        [a,b]=find(ITRI==node1(i));
        j=i;
        [c,b]=find(ITRI(a,:)==node2(j));
        if ~isempty(c)
            jitri=a(c);
            ITRI(ITRI==node1(i))=node2(j);
            ITRI(jitri,:)=[];
            if node1(i)==1
                VER=VER(2:end,:);
            elseif node1(i)==length(VER)
                VER=VER(1:end-1,:);
            else
                VER=[VER(1:node1(i)-1,:);VER(node1(i)+1:end,:)];
            end
            a=find(ITRI(:,1)>node1(i));
            ITRI(a,1)=ITRI(a,1)-1;
            a=find(ITRI(:,2)>node1(i));
            ITRI(a,2)=ITRI(a,2)-1;
            a=find(ITRI(:,3)>node1(i));
            ITRI(a,3)=ITRI(a,3)-1;
        end
    end

else

    for i=1:length(node1)
        [a,b]=find(ITRI==node1(i));
        for j=1:length(node2)
            [c,b]=find(ITRI(a,:)==node2(j));
            if ~isempty(c)
                jitri=a(c);
                ITRI(ITRI==node1(i))=node2(j);
                ITRI(jitri,:)=[];
                if node1(i)==1
                    VER=VER(2:end,:);
                elseif node1(i)==length(VER)
                    VER=VER(1:end-1,:);
                else
                    VER=[VER(1:node1(i)-1,:);VER(node1(i)+1:end,:)];
                end
                a=find(ITRI(:,1)>node1(i));
                ITRI(a,1)=ITRI(a,1)-1;
                a=find(ITRI(:,2)>node1(i));
                ITRI(a,2)=ITRI(a,2)-1;
                a=find(ITRI(:,3)>node1(i));
                ITRI(a,3)=ITRI(a,3)-1;
            end
        end
    end
end
 