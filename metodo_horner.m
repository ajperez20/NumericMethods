function resultado = metodo_horner(coeficientes, x, cifras_significativas)
    %------------------------------------------------------------------------------
    % MÉTODO DE HORNER para evaluar polinomios de forma eficiente.
    %
    % Sintaxis:
    %   resultado = metodo_horner(coeficientes, x)
    %   resultado = metodo_horner(coeficientes, x, cifras_significativas)
    %
    % Parámetros:
    %   - coeficientes: Vector de coeficientes del polinomio ordenados de mayor a menor grado.
    %                   Ejemplo: Para P(x) = 2x^3 - 3x^2 + 4x - 1, usar [2, -3, 4, -1].
    %   - x:            Valor donde se evaluará el polinomio.
    %   - cifras_significativas (Opcional): Número de cifras significativas para el resultado final.
    %                                       Si este parámetro no se proporciona, el resultado se
    %                                       mostrará con la precisión predeterminada de Octave.
    %
    % Salida:
    %   - resultado:    Valor del polinomio evaluado en x.
    %
    %--------------------------------------------------------------------------------
    % Ejemplos de Uso:
    %   Sin especificar cifras significativas:
    %   metodo_horner([1, -7, 8, -0.35], 1.37);
    %
    %   Con 5 cifras significativas:
    %   metodo_horner([1, -7, 8, -0.35], 1.37, 5);
    %
    %--------------------------------------------------------------------------------

    % Inicialización del resultado
    resultado = coeficientes(1);

    % Iteración para aplicar el Método de Horner
    for i = 2:length(coeficientes)
        resultado = resultado * x + coeficientes(i);
    end

    % --- Mostrar el resultado en la consola ---
    fprintf('---------------------------------\n');
    fprintf('Método de Horner:\n');
    fprintf('Polinomio: ');
    % Mostrar el polinomio de forma legible
    grado = length(coeficientes) - 1;
    for i = 1:length(coeficientes)
        if i == 1
            % Manejo del coeficiente principal para evitar 'x^0' si es de grado 0
            if grado == 0
                fprintf('%g', coeficientes(i));
            else
                fprintf('%gx^%d', coeficientes(i), grado);
            end
        else
            current_grado = grado - (i - 1);
            if current_grado == 0 % Término independiente
                if coeficientes(i) >= 0
                    fprintf(' + %g', coeficientes(i));
                else
                    fprintf(' - %g', abs(coeficientes(i)));
                end
            elseif current_grado == 1 % Término lineal (x^1)
                if coeficientes(i) >= 0
                    fprintf(' + %gx', coeficientes(i));
                else
                    fprintf(' - %gx', abs(coeficientes(i)));
                end
            else
                if coeficientes(i) >= 0
                    fprintf(' + %gx^%d', coeficientes(i), current_grado);
                else
                    fprintf(' - %gx^%d', abs(coeficientes(i)), current_grado);
                end
            end
        end
    end
    fprintf('\n');

    % Determinar el formato de salida basado en si se proporcionaron cifras significativas
    fprintf('Evaluado en x = %g:\n', x); % Mostrar x con la precisión predeterminada de Octave

    % `nargin` devuelve el número de argumentos de entrada con los que se llamó a la función.
    if nargin == 3
        % Si se proporcionaron 3 argumentos, significa que se especificaron las cifras significativas
        formato = sprintf('%%.%dg', cifras_significativas);
        fprintf(['Resultado (con %d cifras significativas) = ', formato, '\n'], cifras_significativas, resultado);
    else
        % Si se proporcionaron menos de 3 argumentos, usar la precisión por defecto de Octave
        fprintf('Resultado (precisión predeterminada de Octave) = %g\n', resultado);
    end
    fprintf('---------------------------------\n');
end
