function [raiz, iteraciones] = metodo_biseccion(f, intervalo, tolerancia)
    % Método de la bisección mejorado que muestra intervalos y reemplazos
    %
    % Parámetros:
    %   f - Función a la que se le busca la raíz (debe ser continua)
    %   intervalo - Vector [a, b] que contiene el intervalo inicial
    %   tolerancia - Tolerancia para el criterio de parada
    %
    % Retorna:
    %   raiz - Aproximación de la raíz
    %   iteraciones - Número de iteraciones realizadas

    a = intervalo(1);
    b = intervalo(2);

    % Verificar el teorema de Bolzano
    fprintf('\n=== Verificación inicial ===\n');
    fprintf('f(a) = f(%.6f) = %.6f\n', a, f(a));
    fprintf('f(b) = f(%.6f) = %.6f\n', b, f(b));

    if f(a) * f(b) >= 0
        error('La función no cambia de signo en el intervalo dado (no se cumple el teorema de Bolzano)');
    else
        fprintf('Se cumple f(a)*f(b) < 0 → existe raíz en [%.6f, %.6f]\n\n', a, b);
    end

    iteraciones = 0;
    raiz_anterior = 0;

    fprintf('Iter\t Intervalo usado\t Punto medio\t f(c)\t\t Error\t\t Reemplazo\n');
    fprintf('===============================================================================\n');

    while true
        raiz = (a + b) / 2;
        error_actual = abs(b - a);

        % Determinar qué variable se reemplazará
        if iteraciones > 0
            if f(a) * f(raiz) < 0
                reemplazo = sprintf('b = %.6f', raiz);
            else
                reemplazo = sprintf('a = %.6f', raiz);
            end
        else
            reemplazo = '---';
        end

        fprintf('%d\t [%.9f, %.9f]\t %.9f\t %.9f\t %.9f\t %s\n', ...
                iteraciones, a, b, raiz, f(raiz), error_actual, reemplazo);

        % Criterio de parada
        if f(raiz) == 0 || error_actual < tolerancia
            break;
        end

        % Actualizar intervalo
        if f(a) * f(raiz) < 0
            b = raiz;
        else
            a = raiz;
        end

        raiz_anterior = raiz;
        iteraciones = iteraciones + 1;
    end

    fprintf('\n=== Resultados finales ===\n');
    fprintf('Raíz aproximada: %.8f\n', raiz);
    fprintf('Número de iteraciones: %d\n', iteraciones);
    fprintf('Último intervalo: [%.8f, %.8f]\n', a, b);
    fprintf('f(raíz) = %.8f\n', f(raiz));
    fprintf('Error final: %.8f\n', error_actual);
end
