function [bias, A, U, B1, var1] = unfold3x3SVD(p, lambda, v)
% p should be between 0 and 0.5, lambda should be between 0 and 1
% Function to compute specific matrices for given p and lambda
% Returns bias, A, U, U2, B1, B2 , var1, and var2.

% Define matrix A based on parameter p
A = [1-p, p, 0; p, 1-2*p, p; 0, p, 1-p];

% Define vector R2 and matrix R, and S

S = [ 1/v(1), 0, 0; 0, 1/v(2), 0; 0, 0, 1/v(3)];
R = [-1, 1, 0; 1, -2, 1; 0, 1, -1];
SR = R*S;

% Corrected calculation for U2 
U = (A' * A + lambda * (SR'*SR)) \ A';
% Compute bias matrixes B1
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



%Grafica lambda/bias
%f=@(lambda) unfold3x3SVD(0, lambda, v);
%lambda_values=[0:0.01:1];
%n=length(lambda_values);for i = 1:n;
%sesgo(i)=f(lambda_values(i));
%end;
%disp(sesgo);plot(lambda_values,sesgo);
%xlabel('lambda');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/lambda');



%Grafica p/bias
%f=@(p) unfold3x3SVD(p, 0.1, v);
%p_values=[0:0.01:0.5];
%n=length(p_values);for i = 1:n;
%sesgo(i)=f(p_values(i));
%end;
%disp(sesgo);plot(p_values,sesgo);
%xlabel('p');
%ylabel('bias = e * B * e''');
%title('Gráfica bias/p');
end