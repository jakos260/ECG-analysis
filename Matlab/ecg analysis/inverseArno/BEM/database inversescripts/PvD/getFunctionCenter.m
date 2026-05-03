function lapl(VER,ITRI)

% bool PGLGeometry::getFaceFunctionCenter( const PVector &values, float d, int i,
%                                          P3DVector &p1, P3DVector &p2, P3DVector &pmid ) const

radius = 10;

[ADJ,DIST]=graphdist(ITRI,VER,4);



% bool  found = false;
% s0    = values(ITRI(:,1));
% s1    = values(ITRI(:,2));
% s2    = values(ITRI(:,3));

for i=0:length(VER)
    
    
    
    
    
end


function []=getFunctionAt(VER,ITRI,node,values,d)

%  if ( getFaceFunctionCenter( centerDist, d, i, p1, p2, pmid ) )
%                     {
%                         double angle = acos( (p1 - pmid).dot( (p2 - pmid) ) /
%                                              ( (p1 - pmid).length() * (p2 - pmid).length() ) );
% 
% [nbver,BVER,nbtri,BTRI]=voisin(length(VER),ITRI);
s0    = values(ITRI(:,1));
s1    = values(ITRI(:,2));
s2    = values(ITRI(:,3));


itri= 1:length(ITRI);
remove = find((s0 > d & s1 > d & s2 > d) | (s0 < d & s1 < d & s2 < d));
itri(remove)=[];
I=ITRI;
I(remove,:)=[]
s0    = values(ITRI(itri,1));
s1    = values(ITRI(itri,2));
s2    = values(ITRI(itri,3));

while ~isempty(I)
    
end



if (
    {
        found = true;

        if ( s0 <= d && s1 >= d && s2 >= d )
        {
            VER(I(curItri,1)
            p1 = mGeom.face(i).vertex0() + (mGeom.face(i).vertex1() - mGeom.face(i).vertex0()) * (d - s0) / (s1 - s0);
            p2 = mGeom.face(i).vertex0() + (mGeom.face(i).vertex2() - mGeom.face(i).vertex0()) * (d - s0) / (s2 - s0);

        }
        else if ( s0 >= d && s1 <= d && s2 <= d)
        {
            p1 = mGeom.face(i).vertex0() + (mGeom.face(i).vertex1() - mGeom.face(i).vertex0()) * (s0 - d) / (s0 - s1);
            p2 = mGeom.face(i).vertex0() + (mGeom.face(i).vertex2() - mGeom.face(i).vertex0()) * (s0 - d) / (s0 - s2);
        }
        else if ( s0 <= d && s1 >= d && s2 <= d)
        {
            p1 = mGeom.face(i).vertex1() + (mGeom.face(i).vertex0() - mGeom.face(i).vertex1()) * (s1 - d) / (s1 - s0);
            p2 = mGeom.face(i).vertex1() + (mGeom.face(i).vertex2() - mGeom.face(i).vertex1()) * (s1 - d) / (s1 - s2);
        }
        else if ( s0 >= d && s1 <= d && s2 >= d)
        {
            p1 = mGeom.face(i).vertex1() + (mGeom.face(i).vertex0() - mGeom.face(i).vertex1()) * (d - s1) / (s0 - s1);
            p2 = mGeom.face(i).vertex1() + (mGeom.face(i).vertex2() - mGeom.face(i).vertex1()) * (d - s1) / (s2 - s1);
        }
        else if ( s0 >= d && s1 >= d && s2 <= d)
        {
            p1 = mGeom.face(i).vertex2() + (mGeom.face(i).vertex0()-mGeom.face(i).vertex2()) * (d - s2) / (s0 - s2);
            p2 = mGeom.face(i).vertex2() + (mGeom.face(i).vertex1()-mGeom.face(i).vertex2()) * (d - s2) / (s1 - s2);
        }
        else if ( s0 <= d && s1 <= d && s2 >= d)
        {
            p1 = mGeom.face(i).vertex2() + (mGeom.face(i).vertex0()-mGeom.face(i).vertex2()) * (s2 - d) / (s2 - s0);
            p2 = mGeom.face(i).vertex2() + (mGeom.face(i).vertex1()-mGeom.face(i).vertex2()) * (s2 - d) / (s2 - s1);
        }
        pmid =  p1+(p2-p1)/2 - ( (mGeom.face(i).vertex1() - mGeom.face(i).vertex0()) * (s1 - s0) +
                                 (mGeom.face(i).vertex2() - mGeom.face(i).vertex1()) * (s2 - s1) +
                                 (mGeom.face(i).vertex0() - mGeom.face(i).vertex2()) * (s0 - s2)).normalized() * d; // d should be the distance to the closest maximum or minimum
    }

    return found;
}
