%testing solve function

f = @(x,t) (pi^2 - 1)*exp(-t)*sin(pi*x);
x = @(x) sin(pi*x);
dbc = [0;0];
w = 11; %number of nodes
abc = 8; %timestep 1/abc

answer = solve(w,abc,x,dbc,2,f);

x = linspace(0,1,w);
xn = linspace(0,1,1000);
analytical_solution = exp(-1)*sin(pi*xn);
plot(xn,analytical_solution) % is this even right?
hold on
plot(x,answer(:,abc+1));
title('N = 11, dt = 1/8 BACKWARDS euler')
legend("Actual Solution","galerkin solution")