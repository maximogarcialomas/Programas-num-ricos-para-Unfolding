function [var, A, U, B1,var1] = unfold_var_3x3SVD(p, lambda,v)
% p should be between 0 and 0.5, lambda should be between 0 and 1, 
% v vector stimate 1 x 3.
% Function to compute specific matrices for given p and lambda
% Returns var, A, U, B1, var1

% Define matrix A based on parameter p
A = [1-p, p, 0; p, 1-2*p, p; 0, p, 1-p];


% Define and matrix R, S and SR
S = [ 1/v(1), 0, 0; 0, 1/v(2), 0; 0, 0, 1/v(3)];
R = [-1, 1, 0; 1, -2, 1; 0, 1, -1];
SR = R*S;

% Corrected calculation for U2 
U = (A' * A + lambda * (SR'*SR)) \ A';

% Compute bias matrixes B1 and B2
B1 = U * A - eye(3);

% Compute variances, choose which one to analize
var1 = U' * U;

%Choose the cordenate to analize the variance:
%var = var1(1,1);
var = var1(2,2);
%var=var1(3,3)

%Or choose covariances: i different to j, both between 1 and 3
%var = var1(i,j);

%Grafica lambda/var(/covar)
%f=@(lambda) unfold_var_3x3SVD(0, lambda,v);
%lambda_values=[0:0.01:1];
%n=length(lambda_values);for i = 1:n;
%var_values(i)=f(lambda_values(i));
%end;
%disp(var_values);plot(lambda_values,var_values);
%xlabel('lambda');
%ylabel('variance');
%title('Gráfica variance/lambda');




%Grafica p/var(/covar)
%f=@(p) unfold_var_3x3SVD(p, 0.1,v);
%p_values=[0:0.01:0.5];
%n=length(p_values);for i = 1:n;
%var_values(i)=f(p_values(i));
%end;
%disp(var_values);plot(p_values,var_values);
%xlabel('p');
%ylabel('variance');
%title('Gráfica variance/p');



end