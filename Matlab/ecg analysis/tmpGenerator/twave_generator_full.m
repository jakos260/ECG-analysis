function S = twave_generator_full(t, dep, rep, p)
    % mode = 5; param = [theta(1), 0, theta(2), 0, theta(3)]; % mode 5 - downslope position correction
    % mode = 6; param = [theta(1), 0, theta(2), 0, theta(3)]; % mode 6 - simply smoothing st - downslope transition
    mode = getsMode.Exp_Spline; param = [0.995, 0, p(1), p(2), p(3)]; % mode 7 - spline downslope
    S = gets(t, dep, rep, param, mode);
end