function color=funcolor(CMAP,fun,minfun,maxfun)
dim=size(CMAP);
ncolors=dim(1);
index=round((fun-minfun)/(maxfun-minfun)*ncolors);
index=max(1,index);
index=min(ncolors,index);
color=CMAP(index,:);
