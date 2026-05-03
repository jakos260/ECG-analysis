% nodesw.m
% radio button callback routine of ecgsim
% nosw=get(ui1r3,'Value');
% switches nodes on heart surface in scene1 
nosw=get(ui1r3,'Value');

if nosw ==1,
 if viscene==1, set(hs ,'MarkerEdgeColor','k'); end
 if viscene==7, set(hs ,'MarkerEdgeColor','k'); end
 if viscene==3 & exist('hs3'), set(hs3,'MarkerEdgeColor','k');    end
 else, 
 if viscene==1,  set(hs ,'MarkerEdgeColor','none'); end
 if viscene==7,  set(hs ,'MarkerEdgeColor','none'); end
 if viscene==3 & exist('hs3'), set(hs3,'MarkerEdgeColor','none'); end
end
