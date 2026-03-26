function y = twave_generator(theta, t)
    p = [0.2, 0, theta(1), 0, theta(2)]; % theta(2)
    dep = 0;
    rep = theta(3) * 500;
    S = gets(t, dep, rep, p, 5);
    sMaxIdx = 36;
    y = zeros(1, 500);
    y(1:500-sMaxIdx+1) = S(sMaxIdx:500); % length(y_meas) hardcoded !
end