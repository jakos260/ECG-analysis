% dot product(s): dot=(R1,R2)
%                    =R1(:,1).*R2(:,1)+R1(:,2).*R2(:,2)+R1(:,3).*R2(:,3);
function dot=dots(R1,R2)
	dot=R1(:,1).*R2(:,1)+R1(:,2).*R2(:,2)+R1(:,3).*R2(:,3);
