function [root, iterations, flag] = secantefunction(f, x0, x1, tol_x, tol_f, max_iter)
    % Metodo de la Secante para encontrar raices de una funcion.
    % [root, iterations, flag] = secantefunction(f, x0, x1, tol_x, tol_f, max_iter)
    %
    % Entradas:
    %   f:       La funcion anonima (@(x) ...) de la que se busca la raiz.
    %   x0, x1:  Las dos aproximaciones iniciales distintas.
    %   tol_x:   Tolerancia para la diferencia absoluta entre iteraciones sucesivas (|x_i+1 - x_i| < tol_x).
    %   tol_f:   Tolerancia para el valor absoluto de la funcion en la aproximacion (|f(x_i+1)| < tol_f).
    %   max_iter: Numero maximo de iteraciones permitidas.
    %
    % Salidas:
    %   root:      La aproximacion de la raiz encontrada.
    %   iterations: Una matriz con los datos de cada iteracion [iteracion, x_i, f(x_i)].
    %   flag:      Un indicador de la razon de la terminacion:
    %              1: Convergencia por tolerancia de la funcion (|f(root)| < tol_f).
    %              2: Convergencia por tolerancia del paso (|x_i+1 - x_i| < tol_x).
    %              3: Se alcanzo el numero maximo de iteraciones sin converger dentro de las tolerancias.
    %             -1: Denominador (f(x1) - f(x0)) es cercano a cero, el metodo falla.

    % Inicializacion
    iterations = []; % Matriz para almacenar [iter, x, f(x)]
    iter = 0;
    flag = 0; % Inicializar flag

    % Validar que x0 y x1 sean distintos
    if x0 == x1
        error('Las aproximaciones iniciales x0 y x1 deben ser distintas.');
    end

    % Evaluar la funcion en los puntos iniciales
    f0 = f(x0);
    f1 = f(x1);

    % Añadir la primera iteracion (iter 0) antes de entrar al bucle principal
    % Representamos x1 y f1 como el estado al inicio de la iteracion 0.
    iterations = [iterations; iter, x1, f1];

    % Bucle principal del metodo de la secante
    % Usamos un bucle while true y rompemos cuando se cumpla un criterio
    while true

        % Verificar si la funcion ya esta cerca de cero en la aproximacion actual (x1)
        if abs(f1) < tol_f
            root = x1;
            flag = 1; % Convergencia por tolerancia de la funcion
            break; % Salir del bucle
        end

        % Verificar la diferencia entre iteraciones sucesivas
        % Este criterio se evalua mejor *antes* de calcular la siguiente aproximacion x2,
        % comparando x1 (actual) con x0 (anterior).
        if abs(x1 - x0) < tol_x && iter > 0 % Asegurarse de que ya se hizo al menos una iteracion
             root = x1;
             flag = 2; % Convergencia por tolerancia del paso
             break; % Salir del bucle
        end

        % Verificar el numero maximo de iteraciones
        if iter >= max_iter
            root = x1; % Devuelve la ultima aproximacion encontrada
            flag = 3; % Se alcanzo el maximo de iteraciones
            break; % Salir del bucle
        end

        % Verificar que el denominador no sea cero o muy cercano a cero
        denominator = f1 - f0;
        if abs(denominator) < 1e-12 % Usamos una pequeña tolerancia para evitar division por cero
            root = x1; % Devuelve la ultima aproximacion valida
            flag = -1; % Metodo fallido por denominador cercano a cero
            disp('Advertencia: El denominador (f(x_i) - f(x_i-1)) es cercano a cero.');
            break; % Salir del bucle
        end

        % Calcular la nueva aproximacion (x2)
        x2 = x1 - f1 * (x1 - x0) / denominator;

        % Actualizar los valores para la siguiente iteracion
        x0 = x1; % x0 se convierte en la anterior x1
        f0 = f1; % f0 se convierte en la anterior f1
        x1 = x2; % x1 se convierte en la nueva aproximacion x2

        % Incrementar el contador de iteraciones
        iter = iter + 1;

        % Evaluar la funcion en la nueva aproximacion x1 (que ahora es x2)
        f1 = f(x1);

        % Guardar los datos de la iteracion actual
        iterations = [iterations; iter, x1, f1];

    end % Fin del bucle while

    % --- Mostrar los resultados de las iteraciones con mejor formato ---
    disp(' '); % Linea en blanco para separacion
    disp('Resultados del Metodo de la Secante');
    disp('--------------------------------------');
    fprintf(' Iteracion |      x_i      |      f(x_i)   \n');
    fprintf('--------------------------------------\n');

    % Iterar sobre la matriz de resultados para imprimir cada fila formateada
    for i = 1:size(iterations, 1)
        fprintf(' %8d  | %12.8f  | %12.8f\n', iterations(i, 1), iterations(i, 2), iterations(i, 3));
    end
    fprintf('--------------------------------------\n');

    % Mostrar el motivo de la terminacion
    switch flag
        case 1
            fprintf('Terminacion: Convergencia por tolerancia de la funcion (|f(root)| < %e).\n', tol_f);
        case 2
             % Si se detuvo por tolerancia de paso, la última diferencia fue |x_última - x_anterior|.
             % La diferencia que causo la parada es iterations(end, 2) - iterations(end-1, 2)
             last_x = iterations(end, 2);
             prev_x = iterations(end-1, 2); % Necesitamos al menos 2 iteraciones para este caso
             fprintf('Terminacion: Convergencia por tolerancia del paso (|x_i - x_i-1| < %e).\n', tol_x);
        case 3
            fprintf('Terminacion: Se alcanzo el numero maximo de iteraciones (%d) sin converger.\n', max_iter);
        case -1
            fprintf('Terminacion: El metodo fallo debido a un denominador cercano a cero.\n');
    end
    fprintf('Raiz aproximada encontrada: %12.8f\n', root);
    fprintf('Valor de la funcion en la raiz aproximada: %12.8f\n', f(root)); % Mostrar f(root) final

end % Fin de la funcion
