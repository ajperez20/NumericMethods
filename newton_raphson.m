function [root, iterations, converged] = newton_raphson(f, df, x0, tol, max_iter)
    % Método de Newton-Raphson para encontrar raíces de una función
    %
    % Parámetros de entrada:
    %   f: Función handle de la función a evaluar (ej. @(x) x^2 - 2)
    %   df: Función handle de la derivada de f
    %   x0: Valor inicial (estimación inicial)
    %   tol: Tolerancia para el criterio de parada (opcional, defecto 1e-6)
    %   max_iter: Número máximo de iteraciones (opcional, defecto 100)
    %
    % Parámetros de salida:
    %   root: Aproximación de la raíz encontrada
    %   iterations: Número de iteraciones realizadas
    %   converged: Booleano indicando si el método convergió
    %
    % La función muestra en pantalla los resultados de cada iteración

    % Establecer valores por defecto para parámetros opcionales
    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5
        max_iter = 100;
    end

    % Inicialización de variables
    x = x0;
    iterations = 0;
    converged = false;

    fprintf('Iteración\t x_n\t\t f(x_n)\t\t Error\n');
    fprintf('==================================================\n');

    % Bucle principal del método
    while iterations < max_iter
        fx = f(x);
        dfx = df(x);

        % Mostrar información de la iteración actual
        if iterations == 0
            error_actual = NaN;
        else
            error_actual = abs(x - x_prev);
        end
        fprintf('%5d\t %12.6f\t %12.6f\t %12.6f\n', iterations, x, fx, error_actual);

        % Verificar criterio de convergencia
        if abs(fx) < tol
            converged = true;
            break;
        end

        % Verificar si la derivada es cero (evitar división por cero)
        if abs(dfx) < eps
            fprintf('Derivada cercana a cero. El método puede no converger.\n');
            break;
        end

        % Guardar el valor actual para el cálculo del error en la siguiente iteración
        x_prev = x;

        % Actualizar la estimación usando la fórmula de Newton-Raphson
        x = x - fx / dfx;

        iterations = iterations + 1;
    end

    % Mostrar resultados finales
    fprintf('\nResultado final:\n');
    fprintf('Raíz aproximada: %.8f\n', x);
    fprintf('f(raíz) = %.8f\n', f(x));
    fprintf('Iteraciones realizadas: %d\n', iterations);

    if converged
        fprintf('El método convergió satisfactoriamente.\n');
    else
        fprintf('El método no convergió en el número máximo de iteraciones.\n');
    end

    root = x;
end
