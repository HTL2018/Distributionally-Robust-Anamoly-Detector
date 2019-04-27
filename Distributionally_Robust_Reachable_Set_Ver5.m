% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 17th March, 2019.

clear all; clc; close all;

%% Problem Data

% System Matrices
N  = 100000;
A  = [0.84  0.23
      -0.47 0.12];
B  = [0.07 -0.32
      0.23 0.58];
C  = [1 0
      2 1];
K  = [1.404 -1.042
      1.842 1.008];
L  = [0.0276   0.0448
      -0.01998 -0.0290];
n  = size(A,1);
m  = size(L,1);
res_cov     = [2.086 0.134
               0.134 2.230];
Sigma_w     = [0.045  -0.011               
               -0.011 0.02];
Sigma_v     = 2*eye(n);
mu_noise    = zeros(n,1); 
mu_residual = zeros(n,1);
alarm_rate  = 0.05; 
B_original  = B;
D           = C'*C;
kappa       = res_cov - Sigma_v;
P           = inv(D)*C'*kappa*C*inv(D);

P       = dare((A-L*C)',zeros(n,n),L*Sigma_v*L'+Sigma_w,eye(n));
res_cov = C*P*C'+Sigma_v;

%% Get the distributionally robust residual threshold.
dr_residual_threshold          = 40; % Calculated from previous work 
chi_squared_residual_threshold = ncx2inv(1-alarm_rate,n,0); % 5.99

%% Compute Reachable Sets 

disp('Computing Reachable Sets for Gaussian System Noise & Sensor Noise');
w = mvnrnd(mu_noise,Sigma_w,N)';  
v = mvnrnd(mu_noise,Sigma_v,N)';

x   = zeros(n,N);
e   = zeros(n,N);
r   = zeros(n,N);
z_t = zeros(N);

for j=2:N            
    x(:,j)   = (A + B*K)*x(:,j-1) - B*K*e(:,j-1) + w(:,j-1);        
    e(:,j)   = (A - L*C)*e(:,j-1) - L*v(:,j-1) + w(:,j-1);
    r(:,j-1) = C*e(:,j-1) + v(:,j-1);        
    z_t(:,j-1) = r(:,j-1)'*inv(res_cov)*r(:,j-1);
end

dr_false_alarm         = sum(z_t>dr_residual_threshold)/(N);
chi_square_false_alarm = sum(z_t>chi_squared_residual_threshold)/(N);

%% 
% Plot Gaussian Case - Quadratic Distance Measure in Residual Space

figure;
% Plot the DR threshold
pointA=[-10 15 dr_residual_threshold];
pointB=[10 -15 dr_residual_threshold];
pointC=[10 15 dr_residual_threshold];
normal = cross(pointA-pointB, pointA-pointC); %# Calculate plane normal
syms x y z
P = [x,y,z];
planefunction = dot(normal, P-pointA);
zplane = solve(planefunction, z);
h(1) = fmesh(zplane, [-20, 20, -20, 20],'EdgeColor','blue');
hold on;
% Plot the Chi Squared threshold
pointA=[-10 15 chi_squared_residual_threshold];
pointB=[10 -15 chi_squared_residual_threshold];
pointC=[10 15 chi_squared_residual_threshold];
normal = cross(pointA-pointB, pointA-pointC); %# Calculate plane normal
syms x y z
P = [x,y,z];
planefunction = dot(normal, P-pointA);
zplane = solve(planefunction, z);
h(2) = fmesh(zplane, [-20, 20, -20, 20],'EdgeColor','red');
hold on;
% Plot the quadratic distance measure
h(3) = plot3(r(1,:),r(2,:), z_t, '.k');
grid on
hold on;
xlabel('$r^{x}_{t}$', 'interpreter', 'latex');
ylabel('$r^{y}_{t}$', 'interpreter', 'latex');
zlabel('Quadratic Distance Measure $z_{t}$', 'interpreter', 'latex');
legend(h(1:3),'DR Threshold','Chi-Squared Threshold', '$z_{t}$', 'Interpreter', 'latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 30);
set(gca,'TickLabelInterpreter','latex')
zlim([0 100])
hold off


