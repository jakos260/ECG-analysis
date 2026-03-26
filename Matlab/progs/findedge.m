% file findedge.m
% function edge=findedge(VER,ITRI,node,markedge,edgecol)
% find the nodes forming the loop at the edge of a triangulated patch VER ITRI,
% starting from the node: node.
% In a closed geometry there should be none.
% If a loop is detected, its nodes are listed in a clockwise order
% markedge: any value, but if
% markedge=1:  edge is plotted:  line(VER(edge,1), VER(edge,2), VER(edge,3),'col','w')
% markedge=2:  nodes are marked: text(VER(i,1), VER(i,2), VER(i,3), '*','col','w')
% markedge=3:  plot and mark nodes
% If a plot is required, a call to triplot should preceed !
% see also: buren,voisin,findloop, buurtris
% A. van Oosterom; date: 080719

function edge=findedge(VER,ITRI,node,markedge,edgecol)

if nargin <4, markedge=0;   end
if nargin <5, edgecol='w';  end

ipause=0; % use while debugging

edge=[];
% first, rule out stray nodes
buur=buren(ITRI,node);

if isempty(buur),
    'node is a loner'
    pause
    return,
end

loop=loopnode(ITRI,node);
start=node(1);
if isempty(loop)==0,
    if loop(1)==node(1),
        edge=[edge loop(1)];
        last=loop(end);
        edge=[edge last];
        while last~=start,
            loop=loopnode(ITRI,last);
            last=loop(end);
            if isempty(find(last==edge))==0, start=last; end
            edge=[edge last];
            node=last;
            setnode;
            
            if ipause==1,
                edge
                pause,
            end
            
        end
    end
    
    if markedge==1||markedge==3;,line(VER(edge,1), VER(edge,2), VER(edge,3),'col',edgecol), end
    if markedge>=2, text(VER(edge,1), VER(edge,2), VER(edge,3), '*','color',edgecol), end
    
end
