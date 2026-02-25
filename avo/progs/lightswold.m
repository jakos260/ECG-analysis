% lightsw.m
% radio button callback routine of triplot
% switches licht of patch hs 
lsw=get(ui1r5,'Value');

if lsw==0,
   set(hs ,'Facelighting','none');
else,
    if ~exist('hlight'), lightpos=[1 0 0];end
    hlight=light('position',lightpos,'style','infinite');
    set(hs,'facelighting','phong','Specularstrength',0.2)
end
