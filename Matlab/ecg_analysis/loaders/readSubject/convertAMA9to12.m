%{
ECG12 = [
    VL - VR; VF - VR; VF - VL;
    VR * 1.5; VL * 1.5; VF * 1.5;
    V1; V2; V3; V4; V5; V6;
]
%}
function AMA_12 = convertAMA9to12(AMA_9)
    % Zabezpieczenie przed pustą macierzą lub nieprawidłowym rozmiarem (innym niż standardowe 9 elektrod)
    if isempty(AMA_9) || size(AMA_9, 1) ~= 9
        AMA_12 = AMA_9;
        return;
    end

    % 1. Wyodrębnienie poszczególnych elektrod kończynowych z pierwszych 3 wierszy
    VR = AMA_9(1, :);
    VL = AMA_9(2, :);
    VF = AMA_9(3, :);

    % Wyodrębnienie elektrod przedsercowych (V1-V6) z pozostałych 6 wierszy
    V_chest = AMA_9(4:9, :);

    % Inicjalizacja macierzy wyjściowej 12xN
    AMA_12 = zeros(12, size(AMA_9, 2));

    % 2. Odprowadzenia dwubiegunowe kończynowe (I, II, III)
    AMA_12(1, :) = VL - VR;      % I
    AMA_12(2, :) = VF - VR;      % II
    AMA_12(3, :) = VF - VL;      % III

    % 3. Odprowadzenia jednobiegunowe wzmocnione (aVR, aVL, aVF) ze współczynnikiem 1.5
    AMA_12(4, :) = 1.5 * VR;     % aVR
    AMA_12(5, :) = 1.5 * VL;     % aVL
    AMA_12(6, :) = 1.5 * VF;     % aVF

    % 4. Odprowadzenia przedsercowe (V1-V6) - przepisane bezpośrednio
    AMA_12(7:12, :) = V_chest;
end