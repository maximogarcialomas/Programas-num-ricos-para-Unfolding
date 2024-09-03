function [ var, U, A, B1, var1] = unfoldnxn_var_TU_lineal(p, lambda, K, n)
% p should be between 0 and 0.5, lambda should be between 0 and 1
% v vector stimate 1 x n. y vector to unfold 1 x n.
% Function to compute specific matrices for given p and lambda
% Returns bias, A, U, U2, B1, B2 , var1, and var2.
%Aqui v=yK

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
U = (A' * A + lambda * (R'*R)) \ (A'+lambda*(R'*R)*K');

% Compute bias matrixes B1
B1 = U * A - eye(n);

var1 = U' * U;
%Compute the correlation
D=sqrt(diag(1./diag(var1)));
corr1=D*var1*D;
%Choose the cordenate to analize the variance or the correlation:
% (i equal to j between 1 and n)
var = var1(1,1);
corr = corr1(1,3);
%Or choose covars: (i different to j, both between 1 and n)
%var = var1(i,j);

%Grafica lambda/var(/covar)
%f=@(lambda) unfoldnxn_var_TU_lineal(0,lambda, K, n);
%lambda_values=[0:0.01:1];
%k=length(lambda_values);for i = 1:k;
%var_values(i)=f(lambda_values(i));
%end;
%disp(var_values);plot(lambda_values,var_values);
%xlabel('lambda');
%ylabel('variance');
%title('Gráfica variance/lambda');




%Grafica p/var(/covar)
%f=@(p) unfoldnxn_var_TU_lineal(p, 0.1, K,n);
%p_values=[0:0.01:0.5];
%k=length(p_values);for i = 1:k;
%var_values(i)=f(p_values(i));
%end;
%disp(var_values);plot(p_values,var_values);
%xlabel('p');
%ylabel('variance');
%title('Gráfica variance/p');
end