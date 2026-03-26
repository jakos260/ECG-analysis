function qmarker(pos,color,r)

[v0,i]=make_sphere(3);
v =bsxfun(@plus,v0*r,pos);
qtriplot(v,i)
qtriplot(['color ' color])