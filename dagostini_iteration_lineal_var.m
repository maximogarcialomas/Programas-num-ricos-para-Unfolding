function [corr, K, U, B] = dagostini_iteration_lineal_var(p, lambda0)
    % Esta función realiza la primera iteración de D'Agostini corregida
    % K: matriz de tamaño n x n
    % lambda0: vector inicial de tamaño p
    % El problema es lineal para la para la primera iteracion
    %con v1=(V*K'*W)*y

    n=length(lambda0);
    K = zeros(n);
    
    % Llenar la matriz usando loops
    for i = 1:n
        for j = 1:n
            if i == j
                if i == 1 || i == n
                    K(i, j) = 1 - p;
                else
                    K(i, j) = 1 - 2 * p;
                end
            elseif abs(i - j) == 1
                K(i, j) = p;
            else
                K(i, j) = 0;
            end
        end
    end
    
    %Crear la matriz inicial del desarrollo V
    V = zeros(n);
    for i = 1:n
        V(i,i)=lambda0(i);
    end
    
    %Hacer los sumatorios de los denominadores de W
    t=zeros(1,n);
    for i = 1:n
        t(i)=K(i,:)*lambda0';
    end
    
    %Hacer W
    W = zeros(n);
    for i = 1:n
        W(i,i)=1/t(i);
    end

    %Construir U
    U=(V*K'*W);
    
    % Compute bias matrixes B
    B = U * K - eye(n);

    var1 = U' * U;
%Compute the correlation
D=sqrt(diag(1./diag(var1)));
corr1=D*var1*D;
%Choose the cordenate to analize the variance or the correlation:
% (i equal to j between 1 and n)
var = var1(1,3);
corr = corr1(1,3);
%Or choose covars: (i different to j, both between 1 and n)
%var = var1(i,j);





    %Grafica p/var(/covar)
    %f=@(p) dagostini_iteration_lineal_var(p, lambda0);
    %p_values=[0:0.01:0.5];
    %k=length(p_values);for i = 1:k;
    %var_values(i)=f(p_values(i));
    %end;
    %disp(var_values);plot(p_values,var_values);
    %xlabel('p');
    %ylabel('variance');
    %title('Gráfica variance/p');




end