function [x, Ab] = gauss_elimination(A, b, precision)
    % gauss_elimination Resuelve el sistema de ecuaciones Ax = b usando eliminación de Gauss con pivoteo parcial.
    %
    % Entradas:
    %   A - Matriz de coeficientes (n x n)
    %   b - Vector de términos independientes (n x 1)
    %   precision - Número de decimales para redondear (opcional)
    %
    % Salida:
    %   x - Vector de soluciones (n x 1)
    %   Ab - Matriz aumentada después de la eliminación

    % Obtener el tamaño de la matriz A
    n = size(A, 1);

    % Combinar A y b en una matriz aumentada
    Ab = [A b];

    % Proceso de eliminación
    for k = 1:n-1
        % Pivoteo parcial: encontrar el índice del máximo elemento en la columna k
        [~, max_index] = max(abs(Ab(k:n, k)));
        max_index = max_index + k - 1; % Ajustar el índice

        % Intercambiar filas si es necesario
        if max_index != k
            Ab([k, max_index], :) = Ab([max_index, k], :);
        end

        % Realizar la eliminación
        for i = k+1:n
            % Calcular el factor de eliminación
            factor = Ab(i, k) / Ab(k, k);
            % Restar la fila k multiplicada por el factor de eliminación
            Ab(i, :) = Ab(i, :) - factor * Ab(k, :);
        end
    end

    % Sustitución hacia atrás para encontrar la solución
    x = zeros(n, 1); % Inicializar el vector de soluciones
    for i = n:-1:1
        % Corregir la multiplicación para que sea un producto escalar
        x(i) = (Ab(i, end) - Ab(i, 1:n) * x) / Ab(i, i);
    end

    % Redondear los resultados si se especifica la precisión
    if nargin > 2
        factor = 10^precision;
        x = round(x * factor) / factor;
        Ab = round(Ab * factor) / factor;
    end
end

