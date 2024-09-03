function [bias, A, U, B1, var1] = unfoldnxn(p, lambda,n)
% p should be between 0 and 0.5, lambda should be between 0 and 1 
% n natural
% Function to compute specific matrices for given p and lambda
% Returns A, U, matriz_sesgo, var1

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

% Inverse of matrix A
%A_inv = inv(A);

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



% Corrected calculation for U
U = (A' * A + lambda * R*R) \ A';

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
%f=@(lambda) unfoldnxn(0,lambda,n);
%lambda_values=[0:0.01:1];
%k=length(lambda_values);for i = 1:k;
%sesgo(i)=f(lambda_values(i));
%end;
%disp(sesgo);plot(lambda_values,sesgo);
%xlabel('lambda');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/lambda');




%Grafica p/bias
%f=@(p) unfoldnxn(p,0.1,n);
%p_values=[0:0.01:0.49];
%k=length(p_values);for i = 1:k;
%sesgo(i)=f(p_values(i));
%end;
%disp(sesgo);plot(p_values,sesgo);
%xlabel('p');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/p');

end