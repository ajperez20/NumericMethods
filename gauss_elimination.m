function [x, Ab] = gauss_elimination(A, b, precision)
    % gauss_elimination Resuelve el sistema de ecuaciones Ax = b usando eliminación de Gauss con pivoteo parcial.
    %
    % Entradas:
    %    A - Matriz de coeficientes (n x n)
    %    b - Vector de términos independientes (n x 1)
    %    precision - Número de decimales para redondear (opcional)
    %
    % Salida:
    %    x - Vector de soluciones (n x 1)
    %    Ab - Matriz aumentada después de la eliminación

    fprintf('--- INICIO DEL MÉTODO DE ELIMINACIÓN DE GAUSS ---\n\n');

    % Obtener el tamaño de la matriz A
    n = size(A, 1);

    % Combinar A y b en una matriz aumentada
    Ab = [A b];
    fprintf('**Matriz aumentada inicial:**\n');
    disp(Ab);
    fprintf('\n');

    % Proceso de eliminación
    for k = 1:n-1
        fprintf('**Paso %d: Eliminación en la columna %d**\n', k, k);
        fprintf('Buscando el pivote en la columna %d...\n', k);

        % Pivoteo parcial: encontrar el índice del máximo elemento en la columna k
        [~, max_index] = max(abs(Ab(k:n, k)));
        max_index = max_index + k - 1; % Ajustar el índice

        % Intercambiar filas si es necesario
        if max_index != k
            fprintf('Se intercambia la fila %d con la fila %d para pivoteo.\n', k, max_index);
            Ab([k, max_index], :) = Ab([max_index, k], :);
            disp(Ab);
        else
            fprintf('No es necesario intercambiar filas (el pivote ya está en la posición correcta).\n');
        end
        fprintf('\n');

        % Realizar la eliminación
        for i = k+1:n
            % Calcular el factor de eliminación
            if Ab(k, k) == 0
                error('Error: División por cero. El pivote es cero. El sistema podría no tener solución única.');
            end
            factor = Ab(i, k) / Ab(k, k);
            fprintf('Para hacer cero el elemento A(%d,%d), multiplicamos la fila %d por %.4f y la restamos a la fila %d.\n', i, k, k, factor, i);
            % Restar la fila k multiplicada por el factor de eliminación
            Ab(i, :) = Ab(i, :) - factor * Ab(k, :);
            disp(Ab);
            fprintf('\n');
        end
        fprintf('Matriz después de la eliminación en la columna %d:\n', k);
        disp(Ab);
        fprintf('-----------------------------------------\n');
    end

    fprintf('**Proceso de eliminación hacia adelante completado.**\n');
    fprintf('**Matriz aumentada triangular superior:**\n');
    disp(Ab);
    fprintf('\n');

    % Sustitución hacia atrás para encontrar la solución
    fprintf('**Inicio de la sustitución hacia atrás:**\n');
    x = zeros(n, 1); % Inicializar el vector de soluciones
    for i = n:-1:1
        % Sumar los términos ya conocidos
        sum_terms = 0;
        for j = i+1:n
            sum_terms = sum_terms + Ab(i, j) * x(j);
        end
        % Calcular x(i)
        if Ab(i, i) == 0
            error('Error: División por cero durante la sustitución hacia atrás. El sistema podría no tener solución única.');
        end
        x(i) = (Ab(i, end) - sum_terms) / Ab(i, i);
        fprintf('Calculando x(%d): (%.4f - (%.4f)) / %.4f = %.4f\n', i, Ab(i, end), sum_terms, Ab(i, i), x(i));
    end
    fprintf('\n');

    % Redondear los resultados si se especifica la precisión
    if nargin > 2
        fprintf('**Redondeando resultados a %d decimales...**\n', precision);
        factor = 10^precision;
        x = round(x * factor) / factor;
        Ab = round(Ab * factor) / factor;
        fprintf('Matriz aumentada redondeada:\n');
        disp(Ab);
        fprintf('Vector de soluciones redondeado:\n');
        disp(x);
    end

    fprintf('--- FIN DEL MÉTODO DE ELIMINACIÓN DE GAUSS ---\n\n');
    fprintf('**Resultados finales:**\n');
    fprintf('**Matriz Aumentada Final (después de la eliminación):**\n');
    disp(Ab);
    fprintf('**Vector de Soluciones (x):**\n');
    disp(x);
end
