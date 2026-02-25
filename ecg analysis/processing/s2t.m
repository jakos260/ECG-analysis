function t = s2t(s)
    t = decumsum(s);
end

function y = decumsum(x)
    y = zeros(size(x));
    for i = 2:length(x)
        y(1,i) = x(1,i-1) - x(1,i) - (x(1,1) - x(1,2));
    end
    y(1,1) =  y(1,2);
end