N=100000;

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

As=0.05;
alpha=chi2inv(1-As,2)
  
R_s=[0.045 -0.011;-0.011 0.02];
mean_s=[0;0];
nu=mvnrnd(mean_s,R_s,N)';


R_m=[2 0;0 2];
mean_m=[0;0];
eta=mvnrnd(mean_m,R_m,N)';

P=dare((A-L*C)',zeros(2,2),L*R_m*L'+R_s,eye(2));


x=zeros(2,N);
e=zeros(2,N);
r=zeros(2,N);
z=zeros(1,N);


Sigma=C*P*C'+R_m; 

for i=2:N
    x(:,i)   = (A+B*K)*x(:,i-1)-B*K*e(:,i-1)+nu(:,i-1);
    e(:,i)   = (A-L*C)*e(:,i-1)-L*(eta(:,i-1))+nu(:,i-1);
    r(:,i-1) =  C*e(:,i-1)+eta(:,i-1);
    z_t(:,i-1) =  r(:,i-1)'*inv(Sigma)*r(:,i-1);
end
plot3(r(1,:),r(2,:),z)

chi_squared_false_alarm = sum(z_t>alpha)/N;
dr_false_alarm          = sum(z_t>40)/N;

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
h(3) = plot3(r(1,:),r(2,:),z_t, '.k');
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



