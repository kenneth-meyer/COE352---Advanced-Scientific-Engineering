% HW2, 1.6 #17
% Kenneth Meyer, klm5375
clear all

u = [1;-3;2];
A = [4 1 1; 1 0 2; 1 2 5];

u'*A*u %first output, shows utAu < 0 for some u
eig(A) %second output, shows eigenvalues are not all > 0

% as shown by example, there exists a u, u=(1,-3,2),
% such that utAu < 0, hence the matrix A fails to 
% have utAu > 0. Additionally, A has a negative
% eigenvalue, so it is not a positive definite matrix,
% which also communicates that A fails to have utAu > 0.