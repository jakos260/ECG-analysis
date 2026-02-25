function [yest, G] = rfunc(x, p, needG, ~)
% [yest, G] = rfunc_gen(x, p, needG, ~)
% Wrapper w stylu rfunc do użycia z quamin_em
% x    – wektor czasu (możesz tu podać np. 1:500)
% p    – wektor parametrów, u nas p == theta = [p3, p5, rep]
% needG – 1: licz G (Jacobian), 0: tylko yest
% ~    – placeholder dla funtype, żeby sygnatura była podobna

    % Upewniamy się, że t jest wierszem, bo generator przyjmuje t jako wiersz
    t = x(:).';   % 1 x N

    % Model: wyjście generatora
    y_model = twave_generator(p, t);   % zakładamy, że zwraca wiersz 1 x N
    yest = y_model(:);           % quamin pracuje na wektorze kolumnowym

    % Jeśli Jacobian nie jest potrzebny:
    if nargout < 2 || needG == 0
        G = [];
        return;
    end

    % Liczymy Jacobian numerycznie (centralne różniczkowanie)
    eps = 1e-3; % 1e-4
    npar = numel(p);
    N = numel(yest);
    G = zeros(N, npar);

    for k = 1:npar
        dp = zeros(size(p));
        dp(k) = eps;

        y_plus  = twave_generator(p + dp, t);
        y_minus = twave_generator(p - dp, t);

        y_plus  = y_plus(:);
        y_minus = y_minus(:);

        G(:, k) = (y_plus - y_minus) / (2 * eps);
    end
end
