function [x, iter] = jacobi_gauss_seidel_relative(A, b, x0, tol, max_iter, method)
% jacobi_gauss_seidel_relative Solves Ax=b using Jacobi/GS with relative error criterion
%
%   Inputs:
%       A       - coefficient matrix (n x n)
%       b       - right-hand side vector (n x 1)
%       x0      - initial guess vector (n x 1)
%       tol     - tolerance for the RELATIVE error criterion (||x(k)-x(k-1)||inf / ||x(k)||inf)
%       max_iter- maximum number of iterations to perform
%       method  - string specifying the method: 'jacobi' or 'gauss-seidel'
%
%   Outputs:
%       x       - solution vector (n x 1)
%       iter    - number of iterations performed
%
%   Note: This function uses the relative error criterion ||x(k)-x(k-1)||inf / ||x(k)||inf < tol.
%         The check starts from iteration 1 and requires ||x(k)||inf > 0.

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
        % === Pre-calcular D_inv y L_plus_U para Jacobi ===
        D_inv = diag(1 ./ diagA);
        L_plus_U = L + U;

        % === Cálculo del Radio Espectral para Jacobi ===
        TJ = D_inv * (-L_plus_U);
        spectral_radius = max(abs(eig(TJ)));

    else % gauss-seidel
        % === Pre-calcular D+L y U para Gauss-Seidel ===
        DL = D + L;

        % === Verificar que (D+L) sea invertible para Gauss-Seidel ===
        if rcond(DL) < eps
             error('Matrix (D+L) is singular or nearly singular. Gauss-Seidel not applicable.');
        end

        % === Cálculo del Radio Espectral para Gauss-Seidel ===
        TGS = DL \ (-U);
        spectral_radius = max(abs(eig(TGS)));
    end
    % ---------------------------------------------------------


    % --- Inicio del bucle iterativo ---
    fprintf('\nEjecutando método %s con criterio de error relativo...\n', method);
    fprintf('Radio Espectral (rho(T)): %g\n', spectral_radius);
    fprintf('Tolerancia para error relativo: %g\n', tol);

    fprintf('\nProgreso de las Iteraciones:\n');

    while iter < max_iter
        x_old = x; % Guarda el valor anterior
        iter = iter + 1; % Incrementa el contador

        % === Paso de iteración para el método seleccionado ===
        if strcmp(method,'jacobi')
            x = D_inv * (b - L_plus_U * x_old);
        else % gauss-seidel
            for i=1:n
                sigma = 0;
                for j=1:i-1
                    sigma = sigma + A(i,j)*x(j);
                end
                for j=i+1:n
                    sigma = sigma + A(i,j)*x_old(j);
                end
                x(i) = (b(i) - sigma)/diagA(i);
            end
        end

        % === Calcular y mostrar normas y error relativo en cada iteración ===
        current_norm_inf = norm(x - x_old, inf); % Norma L_inf de la diferencia
        current_x_norm_inf = norm(x, inf);      % Norma L_inf del x actual
        residual = b - A * x; % Vector residual
        residual_norm_L2 = norm(residual, 2); % Norma L_2 del residual

        % Calcular error relativo (solo si iter >= 1 y ||x(k)||inf > 0)
        relative_error = NaN; % Inicializamos como NaN
        if iter >= 1 && current_x_norm_inf > 0
            relative_error = current_norm_inf / current_x_norm_inf;
        end

        % Mostrar resultado de la iteración
        fprintf('  Iteración %d:\n', iter);
        fprintf('    Vector x(k) =\n');
        disp(x);
        fprintf('    ||x(k) - x(k-1)||inf = %g\n', current_norm_inf);
        if isnan(relative_error)
            fprintf('    ||x(k) - x(k-1)||inf / ||x(k)||inf = NaN (||x(k)||inf = 0)\n');
        else
            fprintf('    Error Relativo = %g\n', relative_error); % Mostrar el error relativo
        end
        fprintf('    ||b - A*x(k)||_2     = %g\n', residual_norm_L2);
        fprintf('  --------------------\n');

        % === Verificar el criterio de parada: Error Relativo < tol ===
        % Solo verificamos si el error relativo fue calculable y es menor que la tolerancia
        if ~isnan(relative_error) && relative_error < tol
            break; % Salir del bucle si converge
        end
    end % --- Fin del bucle iterativo ---

    % === Mensaje de resumen final ===
    fprintf('\n--- Resumen Final de Ejecución (%s) ---\n', method);
    fprintf('Radio Espectral (rho(T%s)): %g\n', method, spectral_radius);
    fprintf('Criterio de Parada Utilizado: Error Relativo ||x(k) - x(k-1)||inf / ||x(k)||inf < tol = %g\n', tol);


    % Calcular métricas finales para el resumen
    final_norm_diff_inf = norm(x - x_old, inf); % L_inf diff al final
    final_x_norm_inf = norm(x, inf);          % L_inf de x final
    final_residual_norm_L2 = norm(b - A*x, 2); % L_2 residual final

    final_relative_error = NaN;
    if iter >= 1 && final_x_norm_inf > 0
        final_relative_error = final_norm_diff_inf / final_x_norm_inf;
    end


    % Determinar el estado basado en la última verificación del criterio
    % Si el bucle se rompió Y se calculó un error relativo < tol
    converged_by_criterion = (iter < max_iter) && (~isnan(final_relative_error)) && (final_relative_error < tol);

    if converged_by_criterion
        fprintf('Estado: CONVERGENCIA EXITOSA (según error relativo).\n');
        fprintf('Razón Teórica: rho(T) < 1.\n');
        fprintf('Razón Práctica: El error relativo (%g) fue menor que la tolerancia (%g) en la iteración %d.\n', final_relative_error, tol, iter);
        fprintf('Iteraciones: %d\n', iter);
        fprintf('Norma L2 del Residual Final: %g\n', final_residual_norm_L2);
    else % No convergió por el criterio relativo antes de max_iter
        fprintf('Estado: NO CONVERGIÓ por el criterio de error relativo dentro del máximo de %d iteraciones.\n', max_iter);
        fprintf('Último Error Relativo Calculado: ');
         if isnan(final_relative_error)
            fprintf('NaN (||x(%d)||inf = 0)\n', iter);
        else
            fprintf('%g (>= tol = %g)\n', final_relative_error, tol);
        end
        fprintf('\n'); % Salto de línea para separar

        if spectral_radius >= 1
             fprintf('Razón Teórica: rho(T) >= 1. El método es teóricamente divergente o al borde de la convergencia, lo cual impide o dificulta la convergencia por cualquier criterio.\n');
             fprintf('Norma L2 del Residual Final: %g (puede ser grande si diverge)\n', final_residual_norm_L2);
        else % spectral_radius < 1 pero no convergió por error relativo
             fprintf('Razón Teórica: rho(T) < 1. El método es teóricamente convergente.\n');
             fprintf('Posible Razón Práctica: El número máximo de iteraciones (%d) no fue suficiente para que el error relativo (%g) cayera por debajo de la tolerancia (%g), o la tolerancia es muy estricta para este criterio.\n', max_iter, final_relative_error, tol);
             fprintf('Norma L2 del Residual Final: %g (debería ser pequeña si está convergiendo lentamente)\n', final_residual_norm_L2);
        end

        % Mantenemos la advertencia formal
        warning('JGS:NoConvergence', ...
                'Maximum number of iterations (%d) reached without achieving the relative error tolerance (%g).', ...
                max_iter, tol);
    end
    fprintf('-----------------------------------------------\n');

    % La función retorna x (resultado final) e iter (número total de iteraciones)

end % Fin de la función
