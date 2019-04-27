% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 13th March, 2019.

clear all; clc; close all;

%% Problem Data

% System Matrices
N  = 100;
M  = 200;
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
Sigma_w    = [0.045  -0.011               
              -0.011 0.02];
Sigma_v    = 2*eye(n);
mu_noise   = zeros(n,1);                        
alarm_rate = 0.05; 
B_original = B;

%% Compute Reachable Sets 
noise_param.A           = A;
noise_param.B           = B;
noise_param.mu_noise    = mu_noise;
noise_param.alarm_rate  = alarm_rate;

% Get Worst Case System Noise Distribution
disp('Geting Worst Case System Noise Distribution');
noise_param.Sigma_noise  = Sigma_w;
system_noise_details     = worst_noise_distribution(noise_param);
system_noise_probability = system_noise_details.probabilities;
system_noise_support     = system_noise_details.support;
system_noise_threshold   = system_noise_details.threshold;

% Get Worst Case Sensor Noise Distribution
disp('Getting Worst Case Sensor Noise Distribution');
noise_param.Sigma_noise  = Sigma_v;
sensor_noise_details     = worst_noise_distribution(noise_param);
sensor_noise_probability = sensor_noise_details.probabilities;
sensor_noise_support     = sensor_noise_details.support;
sensor_noise_threshold   = sensor_noise_details.threshold;

% Compute Reachable Sets for Worst Case System Noise & Sensor Noise
disp('Computing Reachable Sets for Worst Case System Noise & Sensor Noise');
for j=1:M
    x_system_noise  = zeros(n,1);
    e_system_noise  = zeros(n,1);
    x_sensor_noise  = zeros(n,1);
    e_sensor_noise  = zeros(n,1);    
    for i=1:N   
        id                      = (j-1)*N + i; 
        w(:,i)                  = system_noise_support(:,sample(system_noise_probability));  
        v(:,i)                  = sensor_noise_support(:,sample(sensor_noise_probability));  
        e_system_noise          = (A - L*C)*e_system_noise + w(:,i);        
        e_sensor_noise          = (A - L*C)*e_sensor_noise - L*v(:,i);
        x_system_noise          = (A + B*K)*x_system_noise - B*K*e_system_noise + w(:,i);
        x_sensor_noise          = (A + B*K)*x_sensor_noise - B*K*e_sensor_noise;        
        x_worst_sys_noise(:,id) = x_system_noise;
        x_worst_sen_noise(:,id) = x_system_noise;
    end
end

% Compute Reachable Sets for Gaussian System Noise & Sensor Noise
disp('Computing Reachable Sets for Gaussian System Noise & Sensor Noise');
for j=1:M        
    x_system_noise  = zeros(n,1);
    e_system_noise  = zeros(n,1);
    x_sensor_noise  = zeros(n,1);
    e_sensor_noise  = zeros(n,1);    
    for i=1:N             
        id                         = (j-1)*N + i; 
        w(:,i)                     = randn(n,1);  
        v(:,i)                     = randn(n,1);  
        e_system_noise             = (A - L*C)*e_system_noise + w(:,i);        
        e_sensor_noise             = (A - L*C)*e_sensor_noise - L*v(:,i);
        x_system_noise             = (A + B*K)*x_system_noise - B*K*e_system_noise + w(:,i);
        x_sensor_noise             = (A + B*K)*x_sensor_noise - B*K*e_sensor_noise;        
        x_gaussian_sys_noise(:,id) = x_system_noise;
        x_gaussian_sen_noise(:,id) = x_system_noise;
    end    
end

% Compute Reachable Sets for Student's T Distributed System Noise & Sensor Noise
disp('Computing Reachable Sets for T Distributed System Noise & Sensor Noise');
for j=1:M    
    x_system_noise  = zeros(n,1);
    e_system_noise  = zeros(n,1);
    x_sensor_noise  = zeros(n,1);
    e_sensor_noise  = zeros(n,1);    
    for i=1:N   
        id                  = (j-1)*N + i;  
        w(:,i)              = rand(n,1); % lognrnd(0,1,n,1);%  trnd(1,n,1); 
        v(:,i)              = rand(n,1); % lognrnd(0,1,n,1);% trnd(1,n,1); 
        e_system_noise      = (A - L*C)*e_system_noise + w(:,i);        
        e_sensor_noise      = (A - L*C)*e_sensor_noise - L*v(:,i);
        x_system_noise      = (A + B*K)*x_system_noise - B*K*e_system_noise + w(:,i);
        x_sensor_noise      = (A + B*K)*x_sensor_noise - B*K*e_sensor_noise;        
        x_t_sys_noise(:,id) = x_system_noise;
        x_t_sen_noise(:,id) = x_system_noise;
    end        
end

% perform the geometric sum between Worst Case System Noise and Sensor Noise reachable sets
disp('Performing the geometric sum between Worst Case System Noise and Sensor Noise reachable sets');
sum_id = id*100;
for i = 1:sum_id
    j    =  floor(rand*id)+1;
    k    =  floor(rand*id)+1;   
    jj   =  floor(rand*id)+1;
    kk   =  floor(rand*id)+1;   
    jjj  =  floor(rand*id)+1;
    kkk  =  floor(rand*id)+1;   
    x_total_worst_noises(:,i)    = x_worst_sys_noise(:,j) + x_worst_sen_noise(:,k);
    x_total_gaussian_noises(:,i) = x_gaussian_sys_noise(:,jj) + x_gaussian_sen_noise(:,kk);
    x_total_t_noises(:,i)        = x_t_sys_noise(:,jjj) + x_t_sen_noise(:,kkk);
end

%% Get the Chi-Squared and DR Ellipsoids

% Obtain the Chi-Squared & Distributionally Robust Ellipsoid for System Noise
noise_param.Sigma_noise = Sigma_w;
noise_param.threshold   = system_noise_threshold;

% Type = 1 - DR Case
disp('Getting DR Ellipsoid for System Noise');
noise_param.type = 1; 
DR_out_param     = noise_bounding_ellipsoid(noise_param);
DR_P_sys_noise   = DR_out_param.P;

% Type = 2 - Chi Squared Case
disp('Getting Chi Squared Ellipsoid for System Noise');
noise_param.type        = 2; 
chi_squared_out_param   = noise_bounding_ellipsoid(noise_param);
Chi_Squared_P_sys_noise = chi_squared_out_param.P;

% Obtain the Chi-Squared & Distributionally Robust Ellipsoid for Sensor Noise
noise_param.Sigma_noise = Sigma_v;
noise_param.threshold   = sensor_noise_threshold;

% Type = 1 - DR Case
disp('Getting DR Ellipsoid for Sensor Noise');
noise_param.type = 1; 
DR_out_param     = noise_bounding_ellipsoid(noise_param);
DR_P_sen_noise   = DR_out_param.P;

% Type = 2 - Chi Squared Case
disp('Getting Chi Squared Ellipsoid for Sensor Noise');
noise_param.type        = 2; 
chi_squared_out_param   = noise_bounding_ellipsoid(noise_param);
Chi_Squared_P_sen_noise = chi_squared_out_param.P;

% Obtain Minkovsky Sum
disp('Computing Minkovsky Sum of Ellipsoids');
dr_total_P          = compute_minkovsky_sum(DR_P_sys_noise,DR_P_sen_noise);
chi_squared_total_P = compute_minkovsky_sum(Chi_Squared_P_sys_noise,Chi_Squared_P_sen_noise);

%% Plotting Code
disp('Plotting The Results');
% Case 1: Worst Case System Noise & Sensor Noise  
% Plot the Chi Squared Ellipsoid
figure;
syms xx yy
ellipse_function  = [xx yy]*chi_squared_total_P*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h(1) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'r','lineWidth',4);
hold on;
% Plot the DR Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*dr_total_P*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h(2) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'b','lineWidth',4);
% Plot the reachable states
h(3) = plot(x_total_worst_noises(1,:),x_total_worst_noises(2,:),'.k');
grid on
hold on;
uistack(h(1), 'top');
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
legend(h(1:3),'Chi-Squared Ellipsoid','DR Ellipsoid', 'Reachable States', 'Interpreter', 'latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 30);
set(gca,'TickLabelInterpreter','latex')
hold off


% Case 2: Gaussian Distributed System Noise & Sensor Noise
% Plot the Chi Squared Ellipsoid
figure;
syms xx yy
ellipse_function  = [xx yy]*chi_squared_total_P*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h(1) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'r','lineWidth',4);
hold on;
% Plot the DR Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*dr_total_P*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h(2) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'b','lineWidth',4);
% Plot the reachable states
h(3) = plot(x_total_gaussian_noises(1,:),x_total_gaussian_noises(2,:),'.k');
grid on
hold on;
uistack(h(1), 'top');
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
legend(h(1:3),'Chi-Squared Ellipsoid','DR Ellipsoid', 'Reachable States', 'Interpreter', 'latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 30);
set(gca,'TickLabelInterpreter','latex')
hold off


% Case 3: T Distributed System Noise & Sensor Noise 
% Plot the Chi Squared Ellipsoid
figure;
syms xx yy
ellipse_function  = [xx yy]*chi_squared_total_P*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h(1) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'r','lineWidth',4);
hold on;
% Plot the DR Ellipsoid
syms xx yy
ellipse_function  = [xx yy]*dr_total_P*[xx;yy]-1;
ellipse_function_handle = matlabFunction(ellipse_function);
h(2) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'b','lineWidth',4);
% Plot the reachable states
h(3) = plot(x_total_t_noises(1,:),x_total_t_noises(2,:),'.k');
grid on
hold on;
uistack(h(1), 'top');
xlabel('$x_1$', 'interpreter', 'latex');
ylabel('$x_2$', 'interpreter', 'latex');
legend(h(1:3),'Chi-Squared Ellipsoid','DR Ellipsoid', 'Reachable States', 'Interpreter', 'latex');
a = findobj(gcf, 'type', 'axes');
h = findobj(gcf, 'type', 'line');
set(h, 'linewidth', 4);
set(a, 'linewidth', 4);
set(a, 'FontSize', 30);
set(gca,'TickLabelInterpreter','latex')
hold off

% % Case 4
% % Plot Gaussian System Noise Reach States and Ellipses
% % Plot the System Noise Chi Squared Ellipsoid
% figure;
% syms xx yy
% ellipse_function  = [xx yy]*Chi_Squared_P_sys_noise*[xx;yy]-1;
% ellipse_function_handle = matlabFunction(ellipse_function);
% h(1) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'r','lineWidth',4);
% hold on;
% % Plot the System Noise DR Ellipsoid
% syms xx yy
% ellipse_function  = [xx yy]*DR_P_sys_noise*[xx;yy]-1;
% ellipse_function_handle = matlabFunction(ellipse_function);
% h(2) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'b','lineWidth',4);
% hold on;
% % Plot the reachable states
% h(3) = plot(x_gaussian_sys_noise(1,:),x_gaussian_sys_noise(2,:),'.k');
% grid on
% hold on;
% uistack(h(1), 'top');
% xlabel('$x_1$', 'interpreter', 'latex');
% ylabel('$x_2$', 'interpreter', 'latex');
% legend(h(1:3),'Chi-Squared Ellipsoid','DR Ellipsoid','Reachable States', 'Interpreter', 'latex');
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 4);
% set(a, 'linewidth', 4);
% set(a, 'FontSize', 30);
% set(gca,'TickLabelInterpreter','latex')
% hold off
% 
% % Case 5
% % Plot Gaussian Sensor Noise Reach States and Ellipses
% % Plot the Sensor Noise Chi Squared Ellipsoid
% figure;
% syms xx yy
% ellipse_function  = [xx yy]*Chi_Squared_P_sen_noise*[xx;yy]-1;
% ellipse_function_handle = matlabFunction(ellipse_function);
% h(1) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'r','lineWidth',4);
% hold on;
% % Plot the Sensor Noise DR Ellipsoid
% syms xx yy
% ellipse_function  = [xx yy]*DR_P_sen_noise*[xx;yy]-1;
% ellipse_function_handle = matlabFunction(ellipse_function);
% h(2) = fimplicit(ellipse_function_handle,[-40 40 -40 40],'b','lineWidth',4);
% hold on;
% % Plot the reachable states
% h(3) = plot(x_gaussian_sen_noise(1,:),x_gaussian_sen_noise(2,:),'.k');
% grid on
% hold on;
% uistack(h(1), 'top');
% xlabel('$x_1$', 'interpreter', 'latex');
% ylabel('$x_2$', 'interpreter', 'latex');
% legend(h(1:3),'Chi-Squared Ellipsoid','DR Ellipsoid','Reachable States', 'Interpreter', 'latex');
% a = findobj(gcf, 'type', 'axes');
% h = findobj(gcf, 'type', 'line');
% set(h, 'linewidth', 4);
% set(a, 'linewidth', 4);
% set(a, 'FontSize', 30);
% set(gca,'TickLabelInterpreter','latex')
% hold off