function [x, iter] = jacobi_gauss_seidel(A, b, x0, tol, max_iter, method)
% jacobi_gauss_seidel Solves the linear system Ax = b using Jacobi or Gauss-Seidel method
%
%   Inputs:
%       A       - coefficient matrix (n x n)
%       b       - right-hand side vector (n x 1)
%       x0      - initial guess vector (n x 1)
%       tol     - tolerance for stopping criteria (for ||x(k)-x(k-1)||inf)
%       max_iter- maximum number of iterations to perform
%       method  - string specifying the method: 'jacobi' or 'gauss-seidel'
%
%   Outputs:
%       x       - solution vector (n x 1)
%       iter    - number of iterations performed
%
%   Example:
%       A = [4 -1 0 -1 0 0;
%            -1 4 -1 0 -1 0;
%             0 -1 4 0 0 -1;
%            -1 0 0 4 -1 0;
%             0 -1 0 -1 4 -1;
%             0 0 -1 0 -1 4];
%       b = [0;5;0;6;-2;6];
%       x0 = zeros(6,1);
%       tol = 1e-6;
%       max_iter = 100;
%       [x_jacobi, iter_jacobi] = jacobi_gauss_seidel(A, b, x0, tol, max_iter, 'jacobi')
%       [x_gs, iter_gs] = jacobi_gauss_seidel(A, b, x0, tol, max_iter, 'gauss-seidel')

    n = length(b);
    x = x0;
    iter = 0;
    spectral_radius = NaN; % Inicializamos el radio espectral

    % === Modificación 1: Verificar el nombre del método ===
    method = lower(method);
    if ~(strcmp(method,'jacobi') || strcmp(method,'gauss-seidel'))
        error('Method must be either ''jacobi'' or ''gauss-seidel''.');
    end

    % === Modificación 2: Verificar elementos diagonales no nulos ===
    diagA = diag(A);
    if any(diagA == 0)
        error('Matrix A must have non-zero diagonal elements for these iterative methods.');
    end

    % === Calcular Matrices de Iteración y Radio Espectral ===
    D = diag(diagA); % Matriz diagonal D
    L = tril(A, -1); % Matriz estrictamente triangular inferior L
    U = triu(A, 1);  % Matriz estrictamente triangular superior U

    if strcmp(method, 'jacobi')
        % === Modificación 3: Pre-calcular D_inv y L_plus_U para Jacobi ===
        D_inv = diag(1 ./ diagA); % Matriz diagonal con inversos de los diagonales
        L_plus_U = L + U; % A sin la diagonal (L+U)

        % === Cálculo del Radio Espectral para Jacobi ===
        TJ = D_inv * (-L_plus_U); % Matriz de iteración de Jacobi TJ = D^-1 * (-(L+U))
        spectral_radius = max(abs(eig(TJ)));

    else % gauss-seidel
        % === Modificación 4: Pre-calcular D+L y U para Gauss-Seidel ===
        DL = D + L; % Matriz triangular inferior D+L

        % === Verificar que (D+L) sea invertible para Gauss-Seidel ===
        if rcond(DL) < eps
             error('Matrix (D+L) is singular or nearly singular. Gauss-Seidel not applicable.');
        end

        % === Cálculo del Radio Espectral para Gauss-Seidel ===
        TGS = DL \ (-U); % Matriz de iteración de Gauss-Seidel TGS = (D+L)^-1 * (-U)
        spectral_radius = max(abs(eig(TGS)));
    end
    % ---------------------------------------------------------


    % --- Inicio del bucle iterativo ---
    fprintf('\nEjecutando método %s...\n', method); % Mensaje al iniciar
    fprintf('Radio Espectral de la matriz de iteración (rho(T)): %g\n', spectral_radius);

    % === Encabezado para la salida iterativa ===
    fprintf('\nProgreso de las Iteraciones:\n');
    % ------------------------------------------

    while iter < max_iter
        x_old = x; % Guarda el valor anterior de x
        iter = iter + 1; % Incrementa el contador de iteraciones

        % === Paso de iteración para el método seleccionado ===
        if strcmp(method,'jacobi')
            % === Modificación 5: Vectorización de Jacobi ===
            x = D_inv * (b - L_plus_U * x_old);
        else % gauss-seidel
            % === Modificación 6: Lógica de Gauss-Seidel in-place ===
            for i=1:n
                sigma = 0;
                for j=1:i-1
                    sigma = sigma + A(i,j)*x(j);
                end
                for j=i+1:n
                    sigma = sigma + A(i,j)*x_old(j);
                end
                % === Modificación 7: Usar diagA precalculado ===
                x(i) = (b(i) - sigma)/diagA(i);
            end
        end

        % === Calcular y mostrar normas en cada iteración (Formato Organizado) ===
        current_norm_inf = norm(x - x_old, inf); % Norma L_inf de la diferencia (criterio de parada)
        residual = b - A * x; % Vector residual
        residual_norm_L2 = norm(residual, 2); % Norma L_2 del residual

        fprintf('  Iteración %d:\n', iter);
        fprintf('    Vector x(k) =\n');
        disp(x); % Muestra el vector x verticalmente, respeta format long
        fprintf('    ||x(k) - x(k-1)||inf = %g\n', current_norm_inf); % Muestra la norma L_inf de la diferencia
        fprintf('    ||b - A*x(k)||_2     = %g\n', residual_norm_L2); % Muestra la norma L_2 del residual
        fprintf('  --------------------\n'); % Separador para claridad entre iteraciones
        % ---------------------------------------------------------------------------

        % Check for convergence using the infinity norm
        % === Modificación 8: Criterio de parada común ===
        if current_norm_inf < tol
            break; % Salir del bucle si converge
        end
    end % --- Fin del bucle iterativo ---

    % === Mensaje de resumen final ampliado ===
    fprintf('\n--- Resumen Final de Ejecución (%s) ---\n', method);
    fprintf('Radio Espectral (rho(T%s)): %g\n', method, spectral_radius);

    final_norm_diff_inf = norm(x - x_old, inf); % Calcula la norma L_inf al final del bucle
    final_residual_norm_L2 = norm(b - A*x, 2); % Calcula la norma L_2 del residual final

    if final_norm_diff_inf < tol && iter <= max_iter
        fprintf('Estado: CONVERGENCIA EXITOSA.\n');
        fprintf('Razón Teórica: rho(T) < 1.\n');
        fprintf('Razón Práctica: El criterio de parada (||x(k) - x(k-1)||inf = %g) fue menor que la tolerancia (%g).\n', final_norm_diff_inf, tol);
        fprintf('Iteraciones: %d\n', iter);
        fprintf('Norma L2 del Residual Final: %g\n', final_residual_norm_L2); % Muestra el residual final
    else
        fprintf('Estado: NO CONVERGIÓ dentro del máximo de %d iteraciones.\n', max_iter);
        fprintf('Último Criterio de Parada (||x(max_iter) - x(max_iter-1)||inf): %g (>= tol = %g)\n', final_norm_diff_inf, tol);

        if spectral_radius >= 1
             fprintf('Razón Teórica: rho(T) >= 1. El método es teóricamente divergente o al borde de la convergencia.\n');
             fprintf('Norma L2 del Residual Final: %g (puede ser grande si diverge)\n', final_residual_norm_L2);
        else % spectral_radius < 1 pero no convergió
             fprintf('Razón Teórica: rho(T) < 1. El método es teóricamente convergente.\n');
             fprintf('Posible Razón Práctica: El número máximo de iteraciones (%d) no fue suficiente o la tolerancia (%g) es demasiado estricta para la tasa de convergencia.\n', max_iter, tol);
             fprintf('Norma L2 del Residual Final: %g (debería ser pequeña si está convergiendo lentamente)\n', final_residual_norm_L2);
        end

        % Mantenemos la advertencia formal de Octave/MATLAB
        warning('JGS:NoConvergence', ...
                'Maximum number of iterations (%d) reached without achieving the desired tolerance (%g).', ...
                max_iter, tol);
    end
    fprintf('-----------------------------------------------\n');
    % ----------------------------------------------------

    % La función retorna x (resultado final) e iter (número total de iteraciones)

end % Fin de la función
