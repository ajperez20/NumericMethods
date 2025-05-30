function [root, iterations, converged] = newton_raphson(f, df, x0, tol, max_iter)
    % Método de Newton-Raphson mejorado con criterios de parada robustos
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
    x_prev = x0;
    error_actual = NaN;

    fprintf('Iteración\t x_n\t\t f(x_n)\t\t df(x_n)\t Error (|x_n - x_{n-1}|)\t |f(x_n)|\n');
    fprintf('==================================================================================\n');

    % Bucle principal del método
    while iterations < max_iter
        fx = f(x);
        dfx = df(x);

        % Verificar criterios de convergencia (ambos deben cumplirse)
        if iterations > 0
            error_actual = abs(x - x_prev);
            if error_actual < tol && abs(fx) < tol
                converged = true;
                break;
            end
        end

        % Mostrar información de la iteración actual
        fprintf('%5d\t %14.9f\t %14.9f\t %14.9f\t\t %14.9f\t\t %14.9f\n', ...
                iterations, x, fx, dfx, error_actual, abs(fx));

        % Verificar si la derivada es cero (evitar división por cero)
        if abs(dfx) < eps
            fprintf('\n¡Advertencia! Derivada cercana a cero (df(x) = %.2e).\n', dfx);
            break;
        end

        % Guardar el valor actual para el cálculo del error en la siguiente iteración
        x_prev = x;

        % Actualizar la estimación usando la fórmula de Newton-Raphson
        x = x - fx / dfx;

        iterations = iterations + 1;
    end

    % Mostrar resultados finales
    fprintf('\n=== RESULTADOS FINALES ===\n');
    fprintf('Raíz aproximada: %.8f\n', x);
    fprintf('f(raíz) = %.8f\n', f(x));
    fprintf('Último error estimado: %.8f\n', error_actual);
    fprintf('Iteraciones realizadas: %d\n', iterations);
    fprintf('Tolerancia especificada: %.1e\n', tol);

    if converged
        fprintf('\nCONVERGENCIA ALCANZADA:\n');
        fprintf('Se cumplieron ambos criterios:\n');
        fprintf('1. |x_n - x_{n-1}| = %.2e < %.1e\n', error_actual, tol);
        fprintf('2. |f(x_n)| = %.2e < %.1e\n', abs(f(x)), tol);
    else
        fprintf('\nADVERTENCIA: El método no convergió completamente:\n');
        if iterations >= max_iter
            fprintf('- Se alcanzó el máximo de iteraciones (%d)\n', max_iter);
        end
        if abs(dfx) < eps
            fprintf('- Derivada muy cercana a cero (posible divergencia)\n');
        end
        fprintf('\nÚltimos valores calculados:\n');
        fprintf('|x_n - x_{n-1}| = %.2e (requerido < %.1e)\n', error_actual, tol);
        fprintf('|f(x_n)| = %.2e (requerido < %.1e)\n', abs(f(x)), tol);
    end

    root = x;
end
