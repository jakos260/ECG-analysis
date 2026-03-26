function [a, b, xm,ym] = linear_regression(x,y,kleur,linew)
% function [a,b, xm,ym] = linear_regression(x,y,kleur,linew);
% compute and draw linear regression line for elements of vectors y and x
% y=ax+b in figure(fig)
% color: kleur; linewidth: xm and ym are mean(x) and mean(y)

% A. van Oosterom; 060822

[ni,nj]=size(x);
if ni<nj, x=x'; end

[ni,nj]=size(y);
if ni<nj, y=y'; end
n=min(size(x,1),size(y,1)); 

%if nargin<3, fig=1; kleur='b'; linew=1; end
if nargin<3, kleur='b'; linew=1; end
if nargin<4, linew=1; end
xm=mean(x);
ym=mean(y);
xmax=max(x);
xmin=min(x);

% rangex=[xmin xmax];
% rangey=[min(y) max(y)];

M=[x-xm ones(n,1)];
ab=pinv(M)*(y-ym);
a=ab(1);
% b=ym+ab(2);

if nargin>2
% figure(fig);
hold on
plot([xmin xmax],[ym-a*(xm-xmin) ym+a*(xmax-xm)],kleur,'linewidth', linew );
end

b=ym-a*xm;

if linew >=2
    disp([num2str(b,4) ' + ' num2str(a,4) ' x'  ])
    if a>0 
        legend([num2str(b,3) ' + ' num2str(a,2) ' x'  ]);
    elseif a>100
        legend([num2str(b,3) ' + ' num2str(a,3) ' x'  ]);
    else
        legend([num2str(b,3) ' - ' num2str(abs(a),3) ' x'  ]);
    end
end
    