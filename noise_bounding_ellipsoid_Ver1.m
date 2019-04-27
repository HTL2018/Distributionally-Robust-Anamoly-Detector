% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 13th March, 2019.

function noise_output_param = noise_bounding_ellipsoid_Ver1(noise_input_param)

    % problem Data
    A           = noise_input_param.A;
    B           = noise_input_param.B;
    C           = noise_input_param.C;
    K           = noise_input_param.K;
    L           = noise_input_param.L;    
    B_original  = B;
    mu_noise    = noise_input_param.mu_noise;
    Sigma_noise = noise_input_param.Sigma_noise;
    alarm_rate  = noise_input_param.alarm_rate;
    type        = noise_input_param.type;
    noise_kind  = noise_input_param.noise_kind;
    sys_cov     = noise_input_param.sys_cov;
    sen_cov     = noise_input_param.sen_cov;
    a_values    = 0.01:0.1:0.99;    
    n           = size(A,1);

    if type == 1           
        noise_threshold = noise_input_param.threshold;
    else
        noise_threshold = ncx2inv(1-alarm_rate,n,0);
    end

    %% Get Ellipsoid using CVX by Solving Folowing SDP
    
    % Get Estimation Error Bounding Ellipsoid
    disp('Obtaining Estimation Error Bounding Ellipsoid');

    % SDP Problem Data
    % A = A, B = I, R = Sigma^{-1}_w / w_bar - System Noise
    % A = A, B = I, R = L*Sigma^{-1}_v*L' / v_bar - Sensor Noise
    
    if noise_kind == 1
        % System Noise
        A_e = A - L*C;
        B_e = eye(n);
        R_e = inv(Sigma_noise)/noise_threshold;
    else
        % Sensor Noise
        A_e = A - L*C;
        B_e = -eye(n);
        R_e = inv(L*Sigma_noise*L')/noise_threshold;
    end    

    % Data Structures
    a_range          = max(size(a_values));
    P_ellipse        = zeros(n,n,a_range);
    Ellipsoid_Volume = zeros(a_range,1);

    for i = 1:a_range

        a = a_values(i)  

        cvx_begin sdp quiet
            variable P(n,n) symmetric       
            minimize (-log_det(P))
            subject to        
                P >= 0; %sdp_tol*eye(n);
                [a*P-A_e'*P*A_e  -A_e'*P*B_e 
                 -B_e'*P*A_e     (1-a)*R_e-B_e'*P*B_e] >= 0;             
        cvx_end

        P_ellipse(:,:,i)  = double(P);
        if(strcmp(cvx_status,'Solved'))
            Ellipsoid_Volume(i) = det(inv(P_ellipse(:,:,i)));
        else
            Ellipsoid_Volume(i) = NaN;
        end
        clear P

    end

    % Extract the minimum volume ellipsoid
    [min_volume,min_index] = nanmin(Ellipsoid_Volume);
    
    
    %% Get State Bounding Ellipsoid
    
    disp('Obtaining State Bounding Ellipsoid');
    % SDP Problem Data
    % p - #of outputs
    % A = A, B = -L*sqrtm(residual_var), R = eye(p)/attack_threshold
    
    if noise_kind == 1
        % System Noise
        A_x  = A + B_original*K;
        B_x  = eye(n);
        P_we = P_ellipse(:,:,min_index); 
        P_1  = noise_threshold*sys_cov;
        P_2  = (B_original*K)*inv(P_we)*(B_original*K)';
        P_3  = compute_minkovsky_sum_Ver1(P_1, P_2); 
        R_x  = P_3;
    else
        % Sensor Noise
        A_x = A + B_original*K;
        B_x = -eye(n);
        P_e = P_ellipse(:,:,min_index);
        R_x = inv((B_original*K)*inv(P_e)*(B_original*K)');
    end
    
    % Data Structures
    a_range          = max(size(a_values));
    P_ellipse        = zeros(n,n,a_range);
    Ellipsoid_Volume = zeros(a_range,1);

    for i = 1:a_range

        a = a_values(i)  

        cvx_begin sdp quiet
            variable P(n,n) symmetric       
            minimize (-log_det(P))
            subject to        
                P >= 0; %sdp_tol*eye(n);
                [a*P-A_x'*P*A_x  -A_x'*P*B_x 
                 -B_x'*P*A_x     (1-a)*R_x-B_x'*P*B_x] >= 0;             
        cvx_end

        P_ellipse(:,:,i)  = double(P);
        if(strcmp(cvx_status,'Solved'))
            Ellipsoid_Volume(i) = det(inv(P_ellipse(:,:,i)));
        else
            Ellipsoid_Volume(i) = NaN;
        end
        clear P

    end

    % Extract the minimum volume ellipsoid
    [min_volume,min_index] = nanmin(Ellipsoid_Volume);
    % Box the output parameter
    noise_output_param.P               = P_ellipse(:,:,min_index);
    noise_output_param.min_volume      = min_volume;

end