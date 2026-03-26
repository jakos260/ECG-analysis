function [ret,alpha] = veronline(p,v1,v2)

ret = v1;
s = ( v2 - v1 );
sq = norm(s);

if ( sq ~= 0 )
    alpha = dot((p - v1), s ) / sq;
    %     alpha = qBound( 0.0, alpha, 1.0 );
    ret = v1 + (s * alpha);
    
    
end