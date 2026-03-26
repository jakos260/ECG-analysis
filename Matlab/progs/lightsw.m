% lightsw.m
% radio button callback routine of triplot
% switches light of patch hs 
lsw=get(ui1r5,'Value');

if lsw==0,
   set(hs ,'Facelighting','none'),
else,
   set(hs,'facelighting','phong','Specularstrength',0.2);
end
