function [theta_opt, y_fit, final_error, iter] = twave_optimize_quamin(y_meas, t)

    % Zakładam, że sygnał ma długość 500 (jak w generatorze)
    if nargin < 2
        t = 1:length(y_meas);
    end

    % Długość na wszelki wypadek przycinamy do 500
    N = min(500, length(y_meas));
    t = t(1:N);
    y_meas = y_meas(1:N);

    % wektor x do LM – użyjemy po prostu indeksów czasowych
    x = t(:);            % kolumna
    y = y_meas(:);       % kolumna

    % theta = [p3, p5, rep]
    theta0 = [-0.1, 0.1, 0.2];   % startowa zgadywanka

    % Marquardt dla generatora
    [theta_opt, iter] = quamin(x, y, theta0, @twave_generator);

    % Dopasowany sygnał
    y_fit = twave_generator(theta_opt, t);

    % Ujednolicenie kształtu
    y_fit = y_fit(:);
    y = y(:);
    M = min(length(y_fit), length(y));

    y_fit = y_fit(1:M);
    y = y(1:M);

    % Błąd dopasowania (norma 2 reszty)
    final_error = rmse(y, y_fit);
end
