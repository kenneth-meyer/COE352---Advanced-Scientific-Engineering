% Trapezoidal rule function for hw3

% output is supposed to be the function; aka the evaluation of u at each pt
function [u] = trapezoidal_rule(dt,f,t_max)
    % note: this function always computes the rule for t = 0 to t=4
    % u'' + u = -cos(pi*t^2)*(4*pi^2*t^2-1) - 2*pi*sin(pi*t^2)
    
    %note: these initial conditions are given. can generalize...
    x = [0;1];
    A = [1,dt/2;-dt/2,1];
    B = [1,-dt/2;dt/2,1];
    q = [0;1]; %given; initial conditions
    c = [dt/2;0];
    
    u = x(2);
    
    %need to make sure i'm not off by one
    for i = dt:dt:t_max
        %not sure how efficient this is; computing inverse every time
        x = A\(B*q + c*(f(i+dt) + f(i)));
        u(length(u)+1) = x(2);
        q = x;
    end 
end