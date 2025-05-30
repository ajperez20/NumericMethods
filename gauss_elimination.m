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

    fprintf('\n==================================================\n');
    fprintf('       MÉTODO DE ELIMINACIÓN DE GAUSS\n');
    fprintf('            CON PIVOTEO PARCIAL\n');
    fprintf('==================================================\n\n');

    % Obtener el tamaño de la matriz A
    n = size(A, 1);

    % Mostrar información inicial
    fprintf('DIMENSIÓN DEL SISTEMA: %d ecuaciones con %d incógnitas\n\n', n, n);
    fprintf('MATRIZ DE COEFICIENTES INICIAL:\n');
    disp(A);
    fprintf('\nVECTOR DE TÉRMINOS INDEPENDIENTES:\n');
    disp(b);
    fprintf('\n');

    % Combinar A y b en una matriz aumentada
    Ab = [A b];
    fprintf('MATRIZ AUMENTADA INICIAL [A|b]:\n');
    disp(Ab);
    fprintf('\n');

    % Proceso de eliminación
    fprintf('--------------------------------------------------\n');
    fprintf('FASE 1: ELIMINACIÓN HACIA ADELANTE CON PIVOTEO\n');
    fprintf('--------------------------------------------------\n\n');

    for k = 1:n-1
        fprintf('ETAPA %d: ELIMINACIÓN EN COLUMNA %d\n', k, k);
        fprintf('--------------------------------\n');

        % Mostrar submatriz actual
        fprintf('Submatriz actual desde la fila %d:\n', k);
        disp(Ab(k:end, k:end));
        fprintf('\n');

        % Pivoteo parcial: encontrar el índice del máximo elemento en la columna k
        [max_val, max_index] = max(abs(Ab(k:n, k)));
        max_index = max_index + k - 1; % Ajustar el índice

        % Mostrar información del pivote
        fprintf('Buscando pivote en columna %d:\n', k);
        fprintf(' - Elemento máximo encontrado: %.4f en fila %d\n', Ab(max_index, k), max_index);
        fprintf(' - Pivote actual: %.4f en fila %d\n', Ab(k, k), k);

        % Intercambiar filas si es necesario
        if max_index != k
            fprintf(' -> Intercambiando fila %d con fila %d\n', k, max_index);
            Ab([k, max_index], :) = Ab([max_index, k], :);
            fprintf('Matriz después del intercambio:\n');
            disp(Ab);
        else
            fprintf(' -> El pivote ya está en la posición correcta\n');
        end
        fprintf('\n');

        % Mostrar el pivote seleccionado
        fprintf('PIVOTE SELECCIONADO: %.4f en posición (%d,%d)\n\n', Ab(k, k), k, k);

        % Realizar la eliminación
        fprintf('OPERACIONES DE ELIMINACIÓN:\n');
        for i = k+1:n
            % Calcular el factor de eliminación
            if Ab(k, k) == 0
                error('ERROR: División por cero. Pivote cero en posición (%d,%d). El sistema podría no tener solución única.', k, k);
            end
            factor = Ab(i, k) / Ab(k, k);

            % Mostrar información de la operación
            fprintf(' - Fila %d = Fila %d - (%.4f) * Fila %d\n', i, i, factor, k);
            fprintf('   Para eliminar el elemento A(%d,%d) = %.4f\n', i, k, Ab(i, k));

            % Restar la fila k multiplicada por el factor de eliminación
            Ab(i, :) = Ab(i, :) - factor * Ab(k, :);

            % Mostrar fila modificada
            fprintf('   Fila %d modificada: ', i);
            fprintf('[ ');
            fprintf('%.4f ', Ab(i, :));
            fprintf(']\n\n');
        end

        fprintf('RESULTADO PARCIAL después de la etapa %d:\n', k);
        disp(Ab);
        fprintf('\n==================================================\n\n');
    end

    fprintf('FASE 1 COMPLETADA: MATRIZ TRIANGULAR SUPERIOR OBTENIDA\n\n');
    fprintf('MATRIZ AUMENTADA FINAL [A|b]:\n');
    disp(Ab);
    fprintf('\n');

    % Sustitución hacia atrás para encontrar la solución
    fprintf('--------------------------------------------------\n');
    fprintf('FASE 2: SUSTITUCIÓN HACIA ATRÁS\n');
    fprintf('--------------------------------------------------\n\n');

    x = zeros(n, 1); % Inicializar el vector de soluciones

    for i = n:-1:1
        % Mostrar la ecuación actual
        fprintf('ECUACIÓN %d: ', i);
        fprintf('%d*x%d ', Ab(i, i), i);
        for j = i+1:n
            fprintf('+ (%d)*x%d ', Ab(i, j), j);
        end
        fprintf('= %d\n', Ab(i, end));

        % Sumar los términos ya conocidos
        sum_terms = 0;
        for j = i+1:n
            sum_terms = sum_terms + Ab(i, j) * x(j);
            fprintf('   + %.4f * %.4f (x%d)\n', Ab(i, j), x(j), j);
        end

        % Calcular x(i)
        if Ab(i, i) == 0
            error('ERROR: División por cero en sustitución hacia atrás. Sistema singular en x(%d).', i);
        end

        x(i) = (Ab(i, end) - sum_terms) / Ab(i, i);

        % Mostrar cálculo detallado
        fprintf('   Cálculo: x%d = (%.4f - %.4f) / %.4f = %.4f\n\n', ...
                i, Ab(i, end), sum_terms, Ab(i, i), x(i));
    end

    % Redondear los resultados si se especifica la precisión
    if nargin > 2
        fprintf('\nREDONDEO DE RESULTADOS a %d decimales...\n', precision);
        factor = 10^precision;
        x_rounded = round(x * factor) / factor;
        Ab_rounded = round(Ab * factor) / factor;

        fprintf('\nSOLUCIONES ORIGINALES:\n');
        disp(x);
        fprintf('\nSOLUCIONES REDONDEADAS:\n');
        disp(x_rounded);

        fprintf('\nMATRIZ AUMENTADA REDONDEADA:\n');
        disp(Ab_rounded);

        % Asignar valores redondeados si se solicita
        x = x_rounded;
        Ab = Ab_rounded;
    end

    fprintf('\n==================================================\n');
    fprintf('            RESULTADOS FINALES\n');
    fprintf('==================================================\n\n');

    fprintf('MATRIZ TRIANGULAR SUPERIOR FINAL [A|b]:\n');
    disp(Ab);
    fprintf('\nSOLUCIÓN DEL SISTEMA (vector x):\n');
    for i = 1:n
        fprintf('   x%d = %.6f\n', i, x(i));
    end

    fprintf('\n==================================================\n');
    fprintf(' FIN DEL MÉTODO DE ELIMINACIÓN DE GAUSS\n');
    fprintf('==================================================\n\n');
end
