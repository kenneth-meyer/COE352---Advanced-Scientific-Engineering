% Kenneth Meyer
% HW3
% 10/18/2020
% Trapezoidal rule problem

% plot the known solution, u(t), compared to v, the approximate sol

dt = 0.025; %varying time step
%function we are using the trapezoidal rule on 
f = @(t) -1*cos(pi*t^2)*(4*pi^2*t^2-1) - 2*pi*sin(pi*t^2);

t = 0:dt:4;
u = cos(pi*t.^2); % QUESTION: "check this"; solve analytically?? or use to compare.

plot(t,u,'g')
hold on

error = zeros(1,4);
n=1;

for dt = [.025, .0125, .00625, .003125]
%for dt = [.1, .01, .001, .0001]
    v = trapezoidal_rule(dt,f,4); %(dt,f(t),t_max)
    t = 0:dt:4;
    plot(t,v)
    hold on
    
    error(n) = abs(cos(pi*4.^2) - v(length(v)));
    n = n+1;
end
legend({'y = cos(pi*t^2)','u(t), dt = .025','u(t), dt = .0125','u(t), dt = .00625','u(t), dt = .003125'},'Location','southwest')
hold off

dt_list = [.025, .0125, .00625, .003125]
error

%this did not work; relationship is clear though. just need to show.k
slope = zeros(2,3); %first row is the change of dt, second row is change of error.
for i=1:3
    slope(1,i) = dt_list(i+1)/dt_list(i);
    slope(2,i) = error(i+1)/error(i);
end
slope

%clearly, when the timestep is halved, the error is quartered. This
%indicates that the error decreases exponentially compared to the decrease
%of dt. This was checked using timesteps that were decreasing by 1/10, and
%the respective errors decreased by 1/100 each time step - hence the error
%decreases at a quadratic rate when compared to dt. (DOUBLE CHECK THIS)

% "how fast is it going to zero compared to dt" numeric or descriptivve
% solution???
