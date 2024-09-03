function [bias, U, A, var1] = unfoldnxnTU_lineal(p, lambda, K, n)
% p should be between 0 and 0.5, lambda should be between 0 and 1
% v vector stimate 1 x n. y vector to unfold 1 x n.
% Function to compute specific matrices for given p and lambda
% Returns bias, A, U, U2, B1, B2 , var1, and var2.
%Aqui v=Ky

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
U = (A' * A + lambda * (R'*R)) \ (A'+lambda*(R'*R)*K);

% Compute bias matrixes B1
B1 = U * A - eye(n);

% Compute variances
var1 = U' * U;

%Define vectors e
% Inicializar el vector
e = zeros(1, n);

% Asignar 1 a la coordenada i
e(3) = 1;
%Calculate bias for e with B1 and B2
bias=e*B1*e';
%bias=e*B2*e';



%Grafica lambda/bias
%f=@(lambda) unfoldnxnTU_lineal(0,lambda,K,n);
%lambda_values=[0:0.01:1];
%k=length(lambda_values);for i = 1:k;
%sesgo(i)=f(lambda_values(i));
%end;
%disp(sesgo);plot(lambda_values,sesgo);
%xlabel('lambda');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/lambda');




%Grafica p/bias
%f=@(p) unfoldnxnTU_lineal(p,0.1, K, n);
%p_values=[0:0.01:0.5];
%k=length(p_values);for i = 1:k;
%sesgo(i)=f(p_values(i));
%end;
%disp(sesgo);plot(p_values,sesgo);
%xlabel('p');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/p');

end