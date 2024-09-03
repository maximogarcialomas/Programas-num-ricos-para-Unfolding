function [x, A] = unfoldnxnTU(p, lambda, v, y,n)
% p should be between 0 and 0.5, lambda should be between 0 and 1
% v vector stimate 1 x n. y vector to unfold 1 x n.
% Function to compute specific matrices for given p and lambda
% Returns bias, A, U, U2, B1, B2 , var1, and var2.

% Inicializar la matriz A
A = zeros(n);
    
    % Llenar la matriz usando loops
    for i = 1:n
        for j = 1:n
            if i == j
                if i == 1 || i == n
                    A(i, j) = 1 - p;
                else
                    A(i, j) = 1 - 2 * p;
                end
            elseif abs(i - j) == 1
                A(i, j) = p;
            else
                A(i, j) = 0;
            end
        end
    end
% Inicializar la matriz R

R = zeros(n);
    
    % Llenar la matriz usando loops
    for i = 1:n
        for j = 1:n
            if i == j
                if i == 1 || i == n
                    R(i, j) = -1;
                else
                    R(i, j) = - 2 ;
                end
            elseif abs(i - j) == 1
                R(i, j) = 1;
            else
                R(i, j) = 0;
            end
        end
    end

% Corrected calculation for x
x = (A' * A + lambda * (R'*R)) \ (A'*y'+lambda*(R'*R)*v');

end