% COE352 HW4
% Kenneth Meyer
% 11/3/2020

% this script is the implementation of a newton method function I defined
% in "newton.m"

guess = [1,1,1,1];
tol = 10^-6;
max_iter = 150;

[sol,e,iter] = newton(guess, tol, max_iter);

fprintf("\nKenneth Meyer hw4 Newton Method\n\n");

fprintf("Initial Guess: [")
fprintf("%g,",guess);
fprintf("]\n");
fprintf("The code converged, and it took %3d iterations\n",iter);

clear all

tol = 10^-6;
max_iter = 150;
guess = [1,-1,1,-1];

[sol,e,iter] = newton(guess, tol, max_iter);
fprintf("Initial Guess: [")
fprintf("%g,",guess);
fprintf("]\n");

fprintf("The code converged and it took %3d iterations\n",iter);

