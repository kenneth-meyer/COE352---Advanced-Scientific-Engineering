%Kenneth Meyer
%klm5375
%12/7/2020
%COE352 Project 2

%testing the code, seeing what M and K are

f = @(x,t) (pi^2 - 1)*exp(-t)*sin(pi*x);

N=11;
n=551;
ts = linspace(0,1,552);

dt = 1/n;
xi = linspace(0,1,N);
    h = xi(2)-xi(1); %gives step size, we assume uniform stepsize
    
    %initialize M and K matricies and F vector (?)
    M = zeros(N,N);
    K = zeros(N,N);
    %f_global = zeros(N-1,n+1);
    f_global = zeros(N,n+1);
    
    %initializing local m and k matricies, size is know due to 
    mlocal = zeros(2,2);
    klocal = zeros(2,2);
    
    %create local to global mapping matrix - first element is 0,1, second
    %element is 1,2, third is 2,3...pattern continues. there are 2 phi
    %functions that exist on each element space.
    local2global_map = [linspace(0,N-2,N-1); linspace(1,N-1,N-1)]';
    
    
    % should this be xi, or eta???im confused
    % parent reference functions (local functions (?))
    phi_1 = @(eta)(1 - eta) / 2; %matrix of these values (?) will be NxN or sum
    phi_2 = @(eta)(1 + eta) / 2;
    %phi = [phi_1,phi_2];
    %derivatives of phi functions; these could be computed using a jacobian
    dphi = [-1/2,1/2]; %dphi1 and dphi2
    deta_dx = 2/h;
    dx_deta = h/2;
    
    % we use the formula below to calculate the quadrature; weights are
    % equal to 1 with the 
    eta_quad_1 = -1/sqrt(3);
    eta_quad_2 = 1/sqrt(3);
    
    %easiest way to set up values for phi1(eta_quad_1)...(2x2 matrix)
    phi_eta = [phi_1(eta_quad_1), phi_1(eta_quad_2);phi_2(eta_quad_1), phi_2(eta_quad_2)];
    % need to check this out...
    
    % get M and K matricies (hardest part)
    % still a little confused on the mapping part/how to do this in
    % code...shouldn't be too bad though
        
    %taken from lecture...
    for k = 1:N-1 %check this, want 11x11 stiffness and mass matricies...
        for l = 1:2
            for m = 1:2
                %mlocal(l,m) = (phi_1(eta_quad_1)*phi_2(eta_quad_1) + phi_1(eta_quad_2)*phi_2(eta_quad_2))*deta_dx;
                mlocal(l,m) = (phi_eta(l,1)*phi_eta(m,1) + phi_eta(l,2)*phi_eta(m,2))*(h);%*dx_deta*2;%*deta_dx;%*dx_deta*2; %added the times 2, idk if it's supposed to be there
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
        end
        
        %eta_quad_1 = eta_quad_1*h;
        %eta_quad_2 = eta_quad_2*h;

        % this is slightly wrong...added negative and it is better tho
        % doesn't get contribution from the RHS, phi2, not sure if this is
        % right or wrong
        %f_global(k,:) = f_global(k,:) - (f(eta_quad_1*h*(k),ts)*phi_1(eta_quad_1) + f(eta_quad_2*h*(k),ts)*phi_1(eta_quad_2))*h/2;%*dx_deta; %added dx_deta
        %f_global(k+1,:) = f_global(k+1,:) - (f(eta_quad_1*h*(k+1),ts)*phi_2(eta_quad_1) + f(eta_quad_1*h*(k+1),ts)*phi_2(eta_quad_2))*h/2;%*dx_deta; %added dx_deta
        %f_global(k,:) = -(f(eta_quad_1*h*(k),ts)*phi_1(eta_quad_1) + f(eta_quad_2*h*(k),ts)*phi_1(eta_quad_2))*h*1.5;
        f_global(k,:) = -(f(eta_quad_1,ts)*phi_1(eta_quad_1) + f(eta_quad_2,ts)*phi_1(eta_quad_2))*(1/8);

        
        %amlost like im forgetting to consider part of it...the *1.5 and
        %negative sign are also complete bs; need to fix.
    end
    %f_global(k+1,:) = -(f(eta_quad_1,ts)*phi_1(eta_quad_1) + f(eta_quad_2,ts)*phi_1(eta_quad_2))*(1/8);

    %f_global(N,:) = 0; % NO IDEA why i put this here...doesn't change much
    %f_global(1,:) = 0;
    f_global
K
M

%do we also somehow get this by applying the bc???think that could be the
%issue....
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

%think I need to apply bc like a did above; before the other steps...

invM = inv(M);
MK = invM*K;

%u(z) is actually a column

u = zeros(N,n+1);

%initial conditions
f_boundary = @(x) sin(pi*x);
u(:,1) = f_boundary(xi);

%dirichlet boundary conditions, idrk how to change
boundaries = [0,0]; %conditions on the boundary
dbc = eye(N);
dbc(1,1) = boundaries(1);
dbc(N,N) = boundaries(2);

size(u);
size(f_global);
size(dbc);
1 - dt*MK; %this was the bane of my existence

%first element is given by initial conditions
for z = 1:n
    t = z*dt;
    u(:,z+1) = u(:,z) - dt*MK*u(:,z) + dt*invM*f_global(:,z); %check f_global
    % now need to apply the boundary condition
    u(:,z+1) = dbc*u(:,z+1);
end
u(:,n+1)
u;
%u(:,n)
%u'
%u(:,1:5)'
%f_global
x = linspace(0,1,11);
xn = linspace(0,1,1000);
sol = exp(-1)*sin(pi*xn);
plot(xn,sol) % is this even right?
hold on
plot(x,u(:,n+1));