function [bv] = unfoldnxnTU_lineal_tradeoff(p, lambda,n,K)
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

%Calculamos el cuadrado del sesgo
B1=B1.^2;

bv=trace(B1+var1);

end
%f=@(lambda) unfoldnxnTU_lineal_tradeoff(p, lambda,n,K);
%lambda_values=[0.1:0.01:1];
%k=length(lambda_values);for i = 1:k;
%bv(i)=f(lambda_values(i));
%end;
%disp(bv);plot(lambda_values,bv);
%xlabel('lambda');
%ylabel('sum(bias^2+var)');
%title('Gráfica sum(bias^2+var)/lambda');