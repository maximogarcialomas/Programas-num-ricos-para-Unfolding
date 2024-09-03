function lambda = dagostini_iteration_nocuad(K ,y, lambda0, num_it)
    % Esta función realiza la iteración de D'Agostini corregida
    % K: matriz de tamaño n x p
    % y: vector de observaciones de tamaño n
    % lambda0: vector inicial de tamaño p
    % num_iterations: número de iteraciones a realizar
    
    [n,p]=size(K);

    % Inicializar lambda con el valor inicial
    lambda = lambda0;
    lambda_new = zeros(p, 1); % Inicializar nuevo lambda
    sum_Kij = zeros(p, 1);
    sum_term = zeros(p, 1);
    sum_Kil_lambda_l = zeros(n, 1);
    % Iterar según el número de iteraciones especificado
    for k=1:num_it
        for i = 1:n
                % Calcular el denominador interno: sum_l (K_il * lambda_l)
                sum_Kil_lambda_l(i) = (K(i, :)* lambda');                
        end
        for j = 1:p
            % Calcular el denominador de la primera fracción: sum_i (K_ij)
            sum_Kij(j) = sum(K(:, j));
            
            % Calcular el segundo término

            
            d=(y./sum_Kil_lambda_l');
            sum_term(j)= d*K(:,j);
            % Actualizar lambda_j
            lambda_new(j) = (lambda(j) / sum_Kij(j)) * sum_term(j);
        end
        
        % Actualizar lambda
        lambda = lambda_new';
    end
end
