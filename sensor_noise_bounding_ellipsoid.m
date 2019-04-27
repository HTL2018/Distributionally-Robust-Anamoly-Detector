% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 13th March, 2019.

function noise_output_param = sensor_noise_bounding_ellipsoid(noise_input_param)

    % problem Data
    A           = noise_input_param.A;
    B           = noise_input_param.B;
    C           = noise_input_param.C;
    K           = noise_input_param.K;
    L           = noise_input_param.L;    
    B_original  = B;    
    Sigma_noise = noise_input_param.Sigma_noise;
    alarm_rate  = noise_input_param.alarm_rate;
    type        = noise_input_param.type;    
    a_values    = 0.01:0.1:0.99;    
    n           = size(A,1);

    if type == 1  
        % DR Case
        noise_threshold = noise_input_param.threshold;
    else
        % Chi-Squared Case
        noise_threshold = ncx2inv(1-alarm_rate,n,0);
    end

    
    %% Get Estimation Error Bounding Ellipsoid
    disp('Obtaining Estimation Error Bounding Ellipsoid');

    % SDP Problem Data
    A_e = A - L*C;
    B_e = eye(n);
    R_e = inv(L*Sigma_noise*L')/noise_threshold;

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
    P_e_v = P_ellipse(:,:,min_index);
    
    %% Get State Bounding Ellipsoid
    
    disp('Obtaining State Bounding Ellipsoid');
    % SDP Problem Data
    A_x = A + B_original*K;
    B_x = eye(n);    
    R_x = inv((B_original*K)*inv(P_e_v)*(B_original*K)');
    
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
    noise_output_param.P = P_ellipse(:,:,min_index);
    
end