% System Matrices
N  = 1000;
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
res_cov     = [2.086 0.134
               0.134 2.230];
Sigma_w     = [0.045  -0.011               
               -0.011 0.02];
Sigma_v     = 2*eye(n);
mu_noise    = zeros(n,1); 
mu_residual = zeros(n,1);
alarm_rate  = 0.05; 
B_original  = B;

D = C'*C;
kappa = res_cov - Sigma_v;
P = inv(D)*C'*kappa*C*inv(D);



dr_residual_threshold          = 40;
chi_squared_residual_threshold = ncx2inv(1-alarm_rate,n,0);

for i=1:N
    w(:,i) = randn(n,1);  
    v(:,i) = randn(n,1);  
end

v_cov = cov(v');
sample_r_cov = C*P*C' + v_cov; 
% Compute Reachable Sets for Worst Case System Noise & Sensor Noise
dr_counter       = 0;
gaussian_counter = 0;
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
        r_worst_residual(:,id)  = C*(e_system_noise+e_sensor_noise) + v(:,i);
        z_worst_measure(id)     = r_worst_residual(:,id)'*inv(res_cov)*r_worst_residual(:,id);
        
        if z_worst_measure(id) > dr_residual_threshold
            dr_counter = dr_counter + 1;
        elseif z_worst_measure(id) > chi_squared_residual_threshold
            gaussian_counter = gaussian_counter + 1;
        end
        
    end
end

worst_case_dr_percent       = dr_counter / id;
worst_case_gaussian_percent = gaussian_counter / id;




% Compute Reachable Sets for Gaussian System Noise & Sensor Noise
disp('Computing Reachable Sets for Gaussian System Noise & Sensor Noise');
dr_counter       = 0;
gaussian_counter = 0;
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
        r_gaussian_residual(:,id)  = C*(e_system_noise+e_sensor_noise) + v(:,i);
        z_gaussian_measure(id)     = r_gaussian_residual(:,id)'*inv(res_cov)*r_gaussian_residual(:,id);
        
        if z_gaussian_measure(id) > dr_residual_threshold
            dr_counter = dr_counter + 1;
        elseif z_gaussian_measure(id) > chi_squared_residual_threshold
            gaussian_counter = gaussian_counter + 1;
        end
        
    end    
end

normal_dr_percent       = dr_counter / id;
normal_gaussian_percent = gaussian_counter / id;
