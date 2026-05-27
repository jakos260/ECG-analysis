function REFLIST=reverse(LIST)
% function REFLIST=reverse(REFLIST)
% outputs rows of LIST in reversed order
[nl jdum]=size(LIST);
if nl==1, LIST=LIST';end
[ni nj]=size(LIST);
LIST=[(ni:-1:1)' LIST];
LIST=sortrows(LIST,1);
REFLIST=LIST(:,2:nj+1);
if nl==1, REFLIST=REFLIST'; end