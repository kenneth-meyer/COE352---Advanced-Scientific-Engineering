
% this solves the heat equation on (x,t) e (0,1)x(0,1)
% N = # nodes
% dt = 1/n (could have framed input differently to make it look nicer)
% ux0 = initial conditions (a boundary condition)
% dbc = [u(0,t), u(1,t)]; initial conditions and dirichlet boundary conditions
% euler_type tells what type of solver we want to use, 1 = fwd, anythingelse = backwards

function [u] = solve(N,n,ux0,bc,euler_type,f)

ts = linspace(0,1,n+1);

dt = 1/n;
xi = linspace(0,1,N);
h = xi(2)-xi(1); %gives step size, uniform step size used for this method
x_eta = @(eta,xi) (eta + 1)*(h/2) + xi;
    
%initialize M and K matricies and F
M = zeros(N,N);
K = zeros(N,N);
f_global = zeros(N,n+1);
    
%initializing local m and k matricies, size is known due to usaage of
%linear basis functions (we are using phi1 and phi2)
mlocal = zeros(2,2);
klocal = zeros(2,2);
    
%local to global mapping matrix, 
local2global_map = [linspace(0,N-2,N-1); linspace(1,N-1,N-1)]';
    
% parent reference functions
phi_1 = @(eta)(1 - eta) / 2;
phi_2 = @(eta)(1 + eta) / 2;
%derivatives of phi functions; these could be computed using a jacobian
dphi = [-1/2,1/2];
deta_dx = 2/h;
dx_deta = h/2;
    
% we use the formula below to calculate the quadrature; weights are
% equal to 1 with the 
eta_quad_1 = -1/sqrt(3);
eta_quad_2 = 1/sqrt(3);
    
%easiest way to set up values for phi1(eta_quad_1)...(2x2 matrix)
phi_eta = [phi_1(eta_quad_1), phi_1(eta_quad_2);phi_2(eta_quad_1), phi_2(eta_quad_2)];
% need to check this out...
    
f_local = zeros(2,n+1);
    
%taken from lecture...
for k = 1:N-1 %check this, want 11x11 stiffness and mass matricies...
    for l = 1:2
        %testing out computing this locally
        y1 = x_eta(eta_quad_1, xi(k+l-1)); %changed x to xi
        y2 = x_eta(eta_quad_2, xi(k+l-1));
        f_local(l,:) = (f(y1,ts)*phi_1(eta_quad_1) + f(y2,ts)*phi_1(eta_quad_2))*h/2;
        for m = 1:2
             mlocal(l,m) = (phi_eta(l,1)*phi_eta(m,1) + phi_eta(l,2)*phi_eta(m,2))*dx_deta;%*dx_deta*2;%*deta_dx;%*dx_deta*2; %added the times 2, idk if it's supposed to be there
             %integral simplifies to below expression
             klocal(l,m) = dphi(l)*deta_dx*dphi(m)*deta_dx*dx_deta*2; %2=evaluated integral from -1 to 1
        end
    end
    for l = 1:2
        global_node = local2global_map(k,l);
        global_node = global_node + 1;
        for m = 1:2
            global_node2 = local2global_map(k,m);
            global_node2 = global_node2 + 1;
            K(global_node,global_node2) = K(global_node,global_node2) + klocal(l,m);
            M(global_node,global_node2) = M(global_node,global_node2) + mlocal(l,m);
        end
        % "el" (spent an hour looking at a mistake, i typed "1" instead of "l" accidentally
        f_global(k+l-1,:) = f_global(k+l-1,:) + f_local(l,:);
    end
        %previous attmepts at f_global. included because i'm still not sure
        %if mine is correct
        %f_global(k,:) = -(f(eta_quad_1,ts)*phi_1(eta_quad_1) + f(eta_quad_2,ts)*phi_1(eta_quad_2))*(1/8);
        %f_global(k,:) = (f(y1,ts)*phi_1(eta_quad_1) + f(y2,ts)*phi_1(eta_quad_2))*h;
end

    %y1 = x_eta(eta_quad_1, x(k+1));
    %y2 = x_eta(eta_quad_2, x(k+1));
    %f_global(k+1,:) = (f(y1,ts)*phi_1(eta_quad_1) + f(y2,ts)*phi_1(eta_quad_2))*h;

    %f_global(k+1,:) = -(f(eta_quad_1,ts)*phi_1(eta_quad_1) + f(eta_quad_2,ts)*phi_1(eta_quad_2))*(1/8);

    
% forcing the M and K matricies to be 0 and 1 on boundary rows/columns
% likely need to alter this slightly to allow for other boundary conditions
M(1,:) = 0;
M(:,1) = 0;
M(N,:) = 0;
M(:,N) = 0;
M(1,1) = 1;
M(N,N) = 1;
M;

%nothing changed with this...
K(1,:) = 0;
K(:,1) = 0;
K(N,:) = 0;
K(:,N) = 0;
K(1,1) = 1;
K(N,N) = 1;
K;

%pre-processing matricies to avoid computing inverse at each step
invM = inv(M);
MK = invM*K;

B = (1/dt)*M + K;
invB = inv(B);

u = zeros(N,n+1);

%initial conditions
%f_boundary = @(x) sin(pi*x);
u(:,1) = ux0(xi);

%dirichlet boundary conditions
boundaries = bc; %defined in a function
dbc = eye(N);
dbc(1,1) = boundaries(1);
dbc(N,N) = boundaries(2);

%forward euler method is type = 1
if euler_type == 1
    %fprintf('aaaa');
    for z = 1:n
        %forward euler
        u(:,z+1) = u(:,z) - dt*MK*u(:,z) + dt*invM*f_global(:,z); %check f_global
        % now need to apply the boundary condition
        u(:,z+1) = dbc*u(:,z+1);
    end
else
    for z = 1:n
        %backward euler
        u(:,z+1) = (1/dt)*invB*M*u(:,z) + invB*f_global(:,z+1); %check f_global, and that it starts at the correct value
        % now need to apply the boundary condition
        u(:,z+1) = dbc*u(:,z+1);
    end
end
end