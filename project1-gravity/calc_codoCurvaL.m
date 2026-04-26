function codo_opt = calc_codoCurvaL(norma_res, norma_sol)
    x = log(norma_res(:));
    y = log(norma_sol(:));

    % extremos
    x1 = x(1); y1 = y(1);
    x2 = x(end); y2 = y(end);

    distancias = zeros(length(x),1);
    for i = 1:length(x)
        % distancia punto-recta
        num = abs((y2 - y1)*x(i) - (x2 - x1)*y(i) + x2*y1 - y2*x1);
        den = sqrt((y2 - y1)^2 + (x2 - x1)^2);
        distancias(i) = num / den;
    end
    
    [~, idx_codo] =max(distancias);
    codo_opt = idx_codo;
end
  
