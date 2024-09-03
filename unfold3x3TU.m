function [x, A] = unfold3x3TU(p, lambda, v, y)
% p should be between 0 and 0.5, lambda should be between 0 and 1
% v vector stimate 1 x 3. y vector to unfold 1 x 3.
% Function to compute specific matrices for given p and lambda
% Returns bias, A, U, U2, B1, B2 , var1, and var2.

% Define matrix A based on parameter p
A = [1-p, p, 0; p, 1-2*p, p; 0, p, 1-p];

% Define and matrix R


R = [-1, 1, 0; 1, -2, 1; 0, 1, -1];

% Corrected calculation for U2 
x = (A' * A + lambda * (R'*R)) \ (A'*y'+lambda*(R'*R)*v');

end