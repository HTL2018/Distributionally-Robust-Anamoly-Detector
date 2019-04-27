% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 13th March, 2019.

function noise_output_param = noise_bounding_ellipsoid(noise_input_param)

    % problem Data
    A           = noise_input_param.A;
    B           = noise_input_param.B;
    mu_noise    = noise_input_param.mu_noise;
    Sigma_noise = noise_input_param.Sigma_noise;
    alarm_rate  = noise_input_param.alarm_rate;
    type        = noise_input_param.type;
    a_values    = 0.01:0.1:0.99;    
    n           = size(A,1);

    if type == 1           
        noise_threshold = noise_input_param.threshold;
    else
        noise_threshold = ncx2inv(1-alarm_rate,n,0);
    end

    %% Get Ellipsoid using CVX by Solving Folowing SDP

    % SDP Problem Data
    % A = A, B = I, R = Sigma^{-1}_w / w_bar - System Noise
    % A = A, B = I, R = Sigma^{-1}_v / v_bar - Sensor Noise

    B = eye(n); 
    R = inv(Sigma_noise)/noise_threshold;

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
                [a*P-A'*P*A  -A'*P*B 
                 -B'*P*A     (1-a)*R-B'*P*B ] >= 0;             
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