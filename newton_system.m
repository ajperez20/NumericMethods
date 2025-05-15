function [x, iter] = newton_system(f, J, x0, tol, max_iter)
    % NEWTON_SYSTEM Método de Newton para sistemas no lineales con impresión de progreso.
    %
    %   [x, iter] = NEWTON_SYSTEM(f, J, x0, tol, max_iter) resuelve un sistema
    %   de ecuaciones no lineales usando el método de Newton, mostrando el progreso.
    %
    %   Parámetros:
    %       f         - Función que toma un vector x y devuelve un vector con las ecuaciones evaluadas.
    %       J         - Función que toma un vector x y devuelve la matriz Jacobiana evaluada en x.
    %       x0        - Vector inicial (punto de inicio).
    %       tol       - Tolerancia para el criterio de convergencia (norma infinita de F(x)).
    %       max_iter  - Número máximo de iteraciones permitidas.
    %
    %   Retorna:
    %       x    - Solución aproximada del sistema.
    %       iter - Número de iteraciones realizadas.
    %
    %   Muestra por pantalla:
    %       - Resultados y error en cada iteración.
    %       - Resultado final y error al converger o al alcanzar el máximo de iteraciones.

    x = x0; % Inicializamos el vector solución
    iter = 0;

    fprintf('--- Iniciando método de Newton ---\n');
    fprintf('Iter |      x      |      Error\n');
    fprintf('----------------------------------\n');

    for iter = 1:max_iter
        F = f(x);           % Evaluar sistema en x
        current_error = norm(F, inf); % Calcular el error (norma infinita de F(x))

        % Imprimir progreso de la iteración
        fprintf('%4d | [%8.4f; %8.4f] | %e\n', iter, x(1), x(2), current_error); % Ajusta el formato si x tiene más dimensiones

        % Verificar convergencia con norma infinita ANTES de calcular delta
        if current_error < tol
            fprintf('----------------------------------\n');
            fprintf('Convergencia alcanzada en %d iteraciones.\n', iter);
            fprintf('Solución final aproximada:\n');
            disp(x);
            fprintf('Error final (norma inf de F(x)): %e\n', current_error);
            return; % Salir de la función
        end

        Jx = J(x);          % Evaluar Jacobiana en x

        % Verificar si la Jacobiana es singular o casi singular
        if rcond(Jx) < eps
             fprintf('----------------------------------\n');
             fprintf('Error: La matriz Jacobiana es singular o casi singular en la iteración %d.\n', iter);
             fprintf('Última solución obtenida:\n');
             disp(x);
             fprintf('Último error (norma inf de F(x)): %e\n', current_error);
            error('Método de Newton falló debido a matriz Jacobiana singular.');
        end

        delta = Jx \ F;     % Resolver sistema lineal J * delta = -F (Octave lo hace como J*delta = F y luego resta delta)
        x = x - delta;      % Actualizar la solución

    end

    % Si el bucle termina sin converger
    fprintf('----------------------------------\n');
    fprintf('Advertencia: Número máximo de %d iteraciones alcanzado.\n', max_iter);
    fprintf('No se alcanzó la convergencia dentro de la tolerancia especificada.\n');
    fprintf('Última solución obtenida:\n');
    disp(x);
    F_final = f(x);
    final_error = norm(F_final, inf);
    fprintf('Último error (norma inf de F(x)): %e\n', final_error);

    % Opcionalmente, puedes lanzar un error si la falta de convergencia es crítica
    % error('No se alcanzó la convergencia después de %d iteraciones.', max_iter);

end
