function [bias, U, A, var1] = unfold3x3TU_lineal(p, lambda, K)
% p should be between 0 and 0.5, lambda should be between 0 and 1
% v vector stimate 1 x 3. y vector to unfold 1 x 3.
% Function to compute specific matrices for given p and lambda
% Returns bias, A, U, U2, B1, B2 , var1, and var2.
%Aqui v=yK

% Define matrix A based on parameter p
A = [1-p, p, 0; p, 1-2*p, p; 0, p, 1-p];

% Define and matrix R


R = [-1, 1, 0; 1, -2, 1; 0, 1, -1];

% Corrected calculation for U2 
U = (A' * A + lambda * (R'*R)) \ (A'+lambda*(R'*R)*K');

% Compute bias matrixes B1 and B2
B1 = U * A - eye(3);
% Compute variances
var1 = U' * U;
%Define vectors e to analize
e=[1, 0, 0];
%e=[0, 1, 0];
%e=[0, 0, 1];
%e=[-1, 0, 1];
%e=[1, -2, 1];

%Calculate bias for e with B1 and B2
bias=e*B1*e';
%bias=e*B2*e';



%Grafica lambda/bias
%f=@(lambda) unfold3x3TU_lineal(0,lambda,K);
%lambda_values=[0:0.01:1];
%n=length(lambda_values);for i = 1:n;
%sesgo(i)=f(lambda_values(i));
%end;
%disp(sesgo);plot(lambda_values,sesgo);
%xlabel('lambda');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/lambda');



%Grafica p/bias
%f=@(p) unfold3x3TU_lineal(p,0.1,K);
%p_values=[0:0.01:0.5];
%n=length(p_values);for i = 1:n;
%sesgo(i)=f(p_values(i));
%end;
%disp(sesgo);plot(p_values,sesgo);
%xlabel('p');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/p');
end