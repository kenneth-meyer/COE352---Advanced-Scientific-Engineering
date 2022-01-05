% COE 352 HW1
% Kenneth Meyer
% klm5375

%% 1.1:3

U = eye(5)-diag(ones(4,1),1);
S = triu(ones(5));

% is this how we have it???
U*S
S
inv(U)

%% 1.1:12

% need to still finish this, but it will be easy.
C=2*eye(4) - diag(ones(3,1),1) - diag(ones(3,1),-1) - diag(ones(1,1),3) - diag(ones(1,1),-3)
[L,U]=lu(C)

%% 1.2:16
fprintf("\n\nSection 1.2:\nProblem 16\n")

% stepsize vectors for k4 and k8
z4 = linspace(0,1,6);
z8 = linspace(0,1,10); %does not really divide evenly

% analytic solution to the function is u(x) = (1/(16pi^2))cos(4pix) -
% (1/(16pi^2))

%k4 matrix
K4 = eye(4)*2 - diag(ones(3,1),1) - diag(ones(3,1),-1);
%k8 matrix
K8 = eye(8)*2 - diag(ones(7,1),1) - diag(ones(7,1),-1);

% K4*u=y ; y is equal to -(1/h^2)cos(4*pi*i*h)
y4 = zeros(1,4);
y8 = zeros(1,8);

h = 0.2;
for i=1:length(y4)
    y4(i) = -(h^2)*cos(4*pi*i*h);
end

h = 1/9;
for i=1:length(y8)
    y8(i) = -(h^2)*cos(4*pi*i*h);
end

y4
y8

%u4 = y4*inv(K4)
%u8 = y8*inv(K8)

u4 = K4*y4'
u8 = K8*y8'

%u4 = K4/y4
%u8 = K8/y8

w = [0];

u4 = [w;u4;w];
u8 = [w;u8;w];
%u4 = [w,u4,w];
%u8 = [w,u8,w];
hold off
fplot(@(x) (1/(16*pi^2))*cos(4*pi*x) - (1/(16*pi^2)),[0,1]);
hold on
plot(z4,u4);
hold on
plot(z8,u8);
legend({'analytical solution','K4','K8'},'Location','Northwest')
title('Analytical and descritized functions for problem 16')
hold off

    
    



%% 1.5:3
fprintf("\n\nSection 1.5:\nProblem 3\n")
% constructing K=K5
K = eye(5)*2 - diag(ones(4,1),1) - diag(ones(4,1),-1);
e = eig(K);
test = 2*ones(5, 1)-2*cos([1:5]*pi/6)';
fprintf("Checking that e agrees with specified line in book:\n");
e-test
fprintf("As the difference between the two is extremely small, e agrees\n");


%% 1.5:16
fprintf("\nProblem 16\n")

A1 = [.6,.4;.4,.6];
A2 = [.6,.9;.1,.6];
A3 = eye(3)*2 - diag(ones(2,1),1) - diag(ones(2,1),-1);

fprintf("A1 eigenvalues:\n");
eig(A1)
fprintf("A2 eigenvalues:\n");
eig(A2)
fprintf("A3 eigenvalues:\n");
eig(A3)





