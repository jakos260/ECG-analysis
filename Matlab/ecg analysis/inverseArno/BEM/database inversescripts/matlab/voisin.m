% voisin.m  (see also: buren, buurtris)
% function [nbver, BVER, nbtri, BTRI, open]=voisin(nver,ITRI)
% for each of the nver nodes of a triangulated
% surface: find 
%               nbver= number of neighouring vertices
%               BVER=  indexes of the neighbouring vertices 
%                      ordered in clock wise fashion 
%                      (NON-CLOSED surfaces permitted)
%               nbtri= number of triangles carrying the node
%               BTRI = indexes of triangles carrying the node
%               open = 1 if local triangulation is
%                      non-closed; else: 0

% comment: 20130307: it seems as if nver in calling specs may be replaced
% by using nver=size(unique(ITRI(:))) in the script; left as is.

function [nbver, BVER, nbtri, BTRI, open]=voisin(nver,ITRI)

dim=size(ITRI);
ntri=dim(1);

nbver=zeros(nver,1);
open =ones(nver,1);
nbtri=zeros(nver,1);

% for each node: find the triangles carrying it
for i=1:ntri,
  for j=1:3,
   ij=ITRI(i,j);
   nbtri(ij)=nbtri(ij)+1;
   BTRI(ij,nbtri(ij))=i;
  end
end

% for each node: find neighboring vertices; arrange
for i=1:nver, loop=1;
% find starting neighbour counter clockwise 
% required for non-closed surface
% set open if required
istart=i;
istop=0;
ind=0;
for idum=1:nbtri(i)+1,loop=2;  
	for j=1:nbtri(i),loop=3;
		for k=1:3,loop=4; 
                  if ITRI(BTRI(i,j),k) == istart,
		    l=k-1; 
		    if l==0, l=3; end
                    ind=ITRI(BTRI(i,j),l);
                    if ind==istop, loop=1; break, 
		    elseif ind ~= i,
                       if istart==i, istop=ind; end 
		       istart=ind; loop=2; break 
                    end
		  end  
		end %loop4               
		if loop<3, break, end
	end %loop3
	if loop<2, break, end
end %loop2
if loop<1, break, end

if loop==1, open(i)=0; end

% start from istart clockwise
nbver(i)=1;
BVER(i,1)=istart;
ingbr=istart;
loop=1;
for idum=1:nbtri(i)+1, loop=2;  
   for j=1:nbtri(i),loop=3;
      for k=1:3,;loop=4;
          if ITRI(BTRI(i,j),k)==ingbr
	     l=icyc(k+1,3);
	     ind=ITRI(BTRI(i,j),l);
	            if ind==BVER(i,1), loop=1; break,
                    elseif ind ~=i
                      nbver(i)=nbver(i)+1;
	              BVER(i,nbver(i))=ind;
	              ingbr=ind;
		      loop=2; break
		    end
          end
      end %loop4		   
   if loop<3, break, end
   end %loop3
if loop<2, break, end
end %loop2
if loop<1, break, end
end %loop1
