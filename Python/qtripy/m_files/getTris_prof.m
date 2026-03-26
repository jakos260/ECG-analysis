function S = getTris_prof(dep, rep, pp, maxt, VER, ITRI )
S = zeros( length(ITRI), maxt );
[p,tr]=triareas(VER,ITRI);
for iFace = 1:length(ITRI)

    v0 = VER(ITRI(iFace,1),:);
    v1 = VER(ITRI(iFace,2),:);
    v2 = VER(ITRI(iFace,3),:);

    tdep0 = dep(ITRI(iFace,1));
    tdep1 = dep(ITRI(iFace,2));
    tdep2 = dep(ITRI(iFace,3));
    triArea = tr(iFace);

    for t = 1:maxt
        if ( t >= tdep0 || t > tdep1 || t >= tdep2 )
            if ( tdep0 <= tdep1 && tdep0 <= tdep2 )
                %s = 1.0 /( 1.0 +  exp(DEPSLOPE  * t-tdep0 ) );
                v01 = (v1-v0) * min(max(0.0, ( t - tdep0 ) / ( tdep1 - tdep0 )), 1.0);
                v02 = (v2-v0) * min(max(0.0, ( t - tdep0 ) / ( tdep2 - tdep0 )), 1.0);

                triangleFraction =norm( cross( v01,v02 ))/ 2.0;
                S(iFace,t) = triangleFraction / triArea;
            elseif tdep1 <= tdep2
                %s = 1.0 /( 1.0 +  exp(DEPSLOPE  * t-tdep1 ) );
                v10 = (v0-v1) * min(max(0.0, ( t - tdep1 ) / ( tdep0 - tdep1 )), 1.0);
                v12 = (v2-v1) * min(max(0.0, ( t - tdep1 ) / ( tdep2 - tdep1 )), 1.0);

                triangleFraction =norm( cross( v10,v12 ))/ 2.0;
                S(iFace,t) = triangleFraction / triArea;
            else
                %s = 1.0; /( 1.0 +  exp(DEPSLOPE  * t-tdep2 ) );
                v20 = (v0-v2) * min(max(0.0, ( t - tdep2 ) / ( tdep0 - tdep2 )),1.0);
                v21 = (v1-v2) * min(max(0.0, ( t - tdep2 ) / ( tdep1 - tdep2 )),1.0);

                triangleFraction =norm( cross(v20,v21 ))/ 2.0;
                S(iFace,t) = triangleFraction / triArea;
            end
            if  S(iFace,t) > 0.999
                trepstart =t;
                S(iFace,t) = 1.0;            
                if ~isempty(rep)
                    trep0 = rep(ITRI(iFace,1));
                    trep1 = rep(ITRI(iFace,2));
                    trep2 = rep(ITRI(iFace,3));
        
                    for t = trepstart:size(S,2)
                    
                        s  = ( 1.0 / ( 1.0 + exp( pp(3) * (t - trep0) ) ) ) * ( 1.0 / ( 1.0 + exp( pp(4) * (t - trep0) ) ) ) +...
                             ( 1.0 / ( 1.0 + exp( pp(3) * (t - trep1) ) ) ) * ( 1.0 / ( 1.0 + exp( pp(4) * (t - trep1) ) ) ) +...
                             ( 1.0 / ( 1.0 + exp( pp(3) * (t - trep2) ) ) ) * ( 1.0 / ( 1.0 + exp( pp(4) * (t - trep2) ) ) );
                       S(iFace,t) = 1-(s / 3.0);         
                    end
                end 
                break
            end
        end
    end
end
