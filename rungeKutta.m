% Assignment 7 problem 6.2.4 - Kenneth Meyer

% u' = -100u + 100sint
f = @(u,t) -100*u + 100*sin(t);

k1 = @(u,t) f(u,t)/2;
k2 = @(u,t,dt) f(u+dt*k1(u,t),t + dt/2)/2;
k3 = @(u,t,dt) f(u+dt*k2(u,t,dt),t + dt/2)/2;
k4 = @(u,t,dt) f(u+2*dt*k3(u,t,dt),t + dt)/2;

g = @(u,t,dt) (1/3)*(k1(u,t) + 2*k2(u,t,dt) + 2*k3(u,t,dt) + k4(u,t,dt));

%testing dt = .0275 and .028 at t = 3 as done in the textbook

dt = [.0275,.028];
for i = 1:2
    t = 0:dt(i):10;
    u = zeros(1,length(t));
    u(1) = 0; %(u(0) = 0)
    
    for j = 1:(length(t)-1)
        u(j+1) = u(j) + g(u(j),t(j),dt(i))*dt(i);
    end
    
    figure(i)
    plot(t,u)
    title(["Runge Kutta, dt = ",num2str(dt(i))]);
    hold off
end



