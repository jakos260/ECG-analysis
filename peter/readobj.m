function shapes=readobj(filename)
fid=fopen(filename);
nshape=1;
vtotal=0;
shape.Name='default';
shape.Vertx=[];
shape.Tri=[];
vstart=0;
while ~feof(fid)
    line=fgetl(fid);
    [cmd,args]=strtok(line);
    switch(cmd)
        case 'o'
            if nshape>0
                shapes{nshape}=shape;
            end
            vstart=vtotal;
            nshape=nshape+1;
            shape.Name=fliplr(deblank(fliplr(args)));
            shape.Vertx=[];
            shape.Tri=[];
        case 'v'
            nV=size(shape.Vertx,1);
            shape.Vertx(nV+1,:)=str2num(args);
            vtotal=vtotal+1;
        case 'f'
            nT=size(shape.Tri,1);
            a=parseslashed(args);
            shape.Tri(nT+1,:)=a(:,1)'-vstart;
    end
end
if nshape>1
    shapes{nshape}=shape;
elseif nshape==1
    shapes=shape;
else
    disp(['no shapes found in file ',filename]);
end


function a=parseslashed(txt)
k=1;
while length(txt)>0
  k2=1;
  [slashed,txt]=strtok(txt);
  while length(slashed)>0
    [item,slashed]=strtok(slashed,'/');
    a(k,k2)=str2num(item);
    k2=k2+1;
  end
  k=k+1;
end

  
