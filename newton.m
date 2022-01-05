% newton method function for hw4

% not entirely sure if we need to generalize or not.
function [sol,e,i] = newton(guess, tol, max_iter)

u = guess';

g_1 = @(u) (u(1) + u(2)^4*u(3));
g_2 = @(u) (-u(1)*u(2) + 2*u(2) + u(3)^3*u(4));
g_3 = @(u) (-2*u(1)^2 + 4*u(3) - u(4));
g_4 = @(u) (2*u(3)^2 + 4*u(4));

%J = jacobian, calculated by hand and using jacobian function in matlab.
J = @(u) ([1, 4*u(2)^3*u(3), u(2)^4,0;-u(2),2 - u(1),3*u(3)^2*u(4),u(3)^3;-4*u(1),0,4,-1;0,0,4*u(3),4]);

i = 0;
e = 1;
while ((i < max_iter) & (e > tol))
    g = [g_1(u);g_2(u);g_3(u);g_4(u)];    
    
    % J(u) is illconditioned (large condition number using cond(J(u))) which leads to 
    % a different number of iterations when using inv(A)*b vs. A\b
    
    %new_u = -inv(J(u))*g + u;
    new_u = -J(u)\g + u;
    
    %calculate error    
    g = [g_1(new_u),g_2(new_u),g_3(new_u),g_4(new_u)];
    e = (g*g')^(1/2); %subtract from zero vector, not needed bc 0 is identity
    
    u = new_u;
    i = i+1;
end
sol = u;
end