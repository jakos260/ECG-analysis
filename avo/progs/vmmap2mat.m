% vmmap2mat.m
% 20070119
% convert Vincent's vmmap  matrix to .mat file

   if ~exist('vmfile') vmfile='vmmap', end
   lchar=size(vmfile,2);
   VM=lsh_read(vmfile);
   VM=VM';
   [nn,nt]=size(VM);
   VM=VM(2:nn,:);
   nn=nn-1;
   fileout=[vmfile '.mat']
   savemat(fileout,VM);


 
