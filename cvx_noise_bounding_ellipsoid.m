% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% MECH 6V29 - Convex Optimization in Systems & Controls
% Project - Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 10th December, 2018.

function noise_output_param = cvx_noise_bounding_ellipsoid(noise_input_param)

    % problem Data
    A          = noise_input_param.A;
    B          = noise_input_param.B;
    mu_noise   = noise_input_param.mu_noise;
    Sigma_w    = noise_input_param.Sigma_w;
    alarm_rate = noise_input_param.alarm_rate;
    a_values   = 0.01:0.1:0.99;
    sdp_tol    = 0.0001;
    n          = size(A,1);

    % Box the input parameters    
    input_param.alarm_rate = alarm_rate;   
    input_param.mu_noise   = mu_noise;   
    input_param.Sigma_w    = Sigma_w;   
    
    % Compute the Noise threshold
    output_param    = compute_noise_threshold(input_param);    
    noise_threshold = output_param.noise_threshold;

    %% Get Ellipsoid using CVX by Solving Folowing SDP

    % SDP Problem Data
    % A = A, B = I, R = Sigma^{-1}_v / noise_threshold

    B = eye(n); 
    R = inv(Sigma_w)/noise_threshold;

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
                 -B'*P*A     (1-a)*R-B'*P*B ] >= 0; %sdp_tol*eye(2*n);                   
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
    P_noise                = P_ellipse(:,:,min_index);

    % Box the output parameter
    noise_output_param.P_noise    = P_noise;
    noise_output_param.min_volume = min_volume;

end