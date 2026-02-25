function S = getTris(dep, rep, maxt, VER, ITRI, plateauSlope, repSlope )
    ntri = size(ITRI,1);
    S = ones(ntri, maxt);
    [~,tr]=triareas(VER,ITRI);

    repSlope = abs(repSlope);
    plateauSlope= abs(plateauSlope)


    for iFace = 1:ntri
        v0 = VER(ITRI(iFace,1),:);
        v1 = VER(ITRI(iFace,2),:);
        v2 = VER(ITRI(iFace,3),:);

        tdep0 = dep(ITRI(iFace,1));
        tdep1 = dep(ITRI(iFace,2));
        tdep2 = dep(ITRI(iFace,3));
        triArea = tr(iFace);

        trepstart = 1;

        for t = 1:maxt
            S(iFace,t) = 0.0;
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
                    
                    S(iFace,t) = 1.0;
                    trepstart = t;
                    break;
                end
            end
        end
        if ~isempty(rep)
            trep0 = rep(ITRI(iFace,1));
            trep1 = rep(ITRI(iFace,2));
            trep2 = rep(ITRI(iFace,3));

            % puoi usare trepstart oppure trepstart+1 come avevi
            for t = (trepstart+1):maxt
                s  = (1.0/(1.0 + exp( plateauSlope*(t - trep0) ))) * (1.0/(1.0 + exp( repSlope*(t - trep0) )));
                s  = s + (1.0/(1.0 + exp( plateauSlope*(t - trep1) ))) * (1.0/(1.0 + exp( repSlope*(t - trep1) )));
                s  = s + (1.0/(1.0 + exp( plateauSlope*(t - trep2) ))) * (1.0/(1.0 + exp( repSlope*(t - trep2) )));
                s  = s / 3.0;

                S(iFace,t) = S(iFace,t) * s;
            end
        end
    end
end


        