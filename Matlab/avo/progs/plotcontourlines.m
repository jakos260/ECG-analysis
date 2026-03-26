% file plotcontourlines.m
% script of triplot_contour
% A. van Oosterom  2014-11-03

% 2011_01_07 stripes option added: overrides the zebra levels
% 2013_11_01 contour option plocked if neither zebra nor stripes are
% specified (do not exist)

tri_nodes_vals = fun(ITRI);
tri_nodes_min  = min(tri_nodes_vals, [], 2);
tri_nodes_max  = max(tri_nodes_vals, [], 2);
absmax         = max(abs([min(VALS) max(VALS)]));

if ~exist('zebra'),
    levels  = linspace(-absmax,absmax,21);
    levels(1) = [];
    levels(end) = [];
elseif zebra >= 0,
    levels = linspace(-absmax,absmax,(zebra+1)*2+1);
    levels(1) = [];
    levels(end) = [];
elseif zebra < 0,
    nlagen=floor(absmax/(-zebra));
    levels=(-nlagen:nlagen)*(-zebra);
    levels(1) = [];
end

if exist('stripes'),
    if ~isempty(stripes),
        levels=stripes;
    end
end

nlevs=size(levels,2);

if ~isempty(nlevs),
    if nlevs>=1,
        
        for ilev=1:nlevs,
            level = levels(ilev);
            use = level>=tri_nodes_min & level<=tri_nodes_max;
            intersect1 = [];
            intersect2 = [];
            
            fun1=[];
            fun2=[];
            
            for itri=find(use)',
                pos  = VER(ITRI(itri,:), :);
                v(1) = tri_nodes_vals(itri,1);
                v(2) = tri_nodes_vals(itri,2);
                v(3) = tri_nodes_vals(itri,3);
                
                la(1) = (level-v(1)) / (v(2)-v(1)+ eps); % abcissa between vertex 1 and 2
                la(2) = (level-v(2)) / (v(3)-v(2)+ eps); % abcissa between vertex 2 and 3
                la(3) = (level-v(3)) / (v(1)-v(3)+ eps); % abcissa between vertex 1 and 2
                
                sel = find(la>=0 & la<=1);        
                
                if size(sel,2)>1;
                    abc(1,:) = pos(1,:) + la(1) * (pos(2,:) - pos(1,:));
                    abc(2,:) = pos(2,:) + la(2) * (pos(3,:) - pos(2,:));
                    abc(3,:) = pos(3,:) + la(3) * (pos(1,:) - pos(3,:));
                    abc=abc(sel,:);
                    if exist('fun0')                   
                       fff(1,:) = fun0(ITRI(itri,1)) + la(1) * (fun0(ITRI(itri,2)) - fun0(ITRI(itri,1)));
                       fff(2,:) = fun0(ITRI(itri,2)) + la(2) * (fun0(ITRI(itri,3)) - fun0(ITRI(itri,2)));
                       fff(3,:) = fun0(ITRI(itri,3)) + la(3) * (fun0(ITRI(itri,1)) - fun0(ITRI(itri,3)));
                       fff=fff(sel,:);
                    end
                    
                    if size(abc,1)>1,
                        if norm3d(abc(1,:)-abc(2,:))< sqrt(eps),
                            abc(2,:)=[];
                            if exist('fun0'), fff(2)=[];end
                            if  size(abc,1)>1 & norm3d(abc(1,:)-abc(2,:))< sqrt(eps),
                                abc(2,:)=[];
                                if exist('fun0'), fff(2)=[];end
                            end
                        end
                    end
                    
                    if size(abc,1)>=2,
                        intersect1 = [intersect1; abc(1,:)];
                        intersect2 = [intersect2; abc(2,:)];
                        if exist('fun0'),
                           fun1 = [fun1; fff(1)];
                           fun2 = [fun2; fff(2)]; 
                        end
                    end
                end
                
            end
            
            % store the details for external reference
          
            cntour(ilev).level = level;
            cntour(ilev).n     = size(intersect1,1);
            cntour(ilev).intersect1 = intersect1;
            cntour(ilev).intersect2 = intersect2;
        end
         
        % collect all different contourlevels and plot them
        intersect1 = [];
        intersect2 = [];
        cntlevel   = [];
        
        for ilev=1:nlevs,
            intersect1 = [intersect1; cntour(ilev).intersect1];
            intersect2 = [intersect2; cntour(ilev).intersect2];
            cntlevel   = [cntlevel; ones(cntour(ilev).n,1) * levels(ilev)];
        end
        
        if isempty(intersect1)==0 & isempty(intersect1)==0 & isempty(cntlevel)==0,
            
            X = [intersect1(:,1) intersect2(:,1)]';
            Y = [intersect1(:,2) intersect2(:,2)]';
            C = [cntlevel(:) cntlevel(:)]'; 
            
            if size(VER,2)>2
                Z = [intersect1(:,3) intersect2(:,3)]';
            else
                Z = zeros(2, length(cntlevel));
            end
            
            hc = [];
            for i=1:max(size(cntlevel)),
                h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
                     'ZData', Z(:,i),'edgecolor',contourcolor);
                hc = [hc; h1];
            end
        end
    end
end
