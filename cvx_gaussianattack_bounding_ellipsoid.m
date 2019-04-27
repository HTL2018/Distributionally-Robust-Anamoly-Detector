% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% MECH 6V29 - Convex Optimization in Systems & Controls
% Distributionally Robust Ellipsoid for Reachable Set Under Zero Alarm Attack
% Date: 9th December, 2018.

function attack_output_param = cvx_gaussianattack_bounding_ellipsoid(attack_input_param)

    % problem Data
    A            = attack_input_param.A;
    B            = attack_input_param.B;
    L            = attack_input_param.L;
    G            = attack_input_param.G;
    B_original   = B;
    mu_attack    = attack_input_param.mu_attack;
    alarm_rate   = attack_input_param.alarm_rate;
    attack_cov   = attack_input_param.attack_cov;
    residual_var = attack_input_param.residual_var;
    a_values     = 0.01:0.1:0.99;
    sdp_tol      = 0.0001;
    n            = size(A,1);
    

    % Box the input parameters    
    input_param.alarm_rate   = alarm_rate;   
    input_param.mu_attack    = mu_attack;   
    input_param.attack_cov   = attack_cov;
    input_param.residual_var = residual_var;  
    
    attack_threshold = ncx2inv(1-alarm_rate,n,0);
    % attack_threshold = attack_threshold + 160;
    
    %% Get Ellipsoid for Zero Alarm Attack - Use CVX to Solve Folowing SDP

    % SDP Problem Data
    % p - #of outputs
    % A = A, B = -L*sqrtm(residual_var), R = eye(p)/attack_threshold
    
    B = -L*sqrtm(residual_var);
    R = eye(n)/attack_threshold;
    
    % Data Structures
    a_range          = max(size(a_values));
    P_ellipse        = zeros(n,n,a_range);
    Ellipsoid_Volume = zeros(a_range,1);

    for i = 1:a_range

        clear p
        a = a_values(i)

        % Solve SDP
        cvx_begin sdp quiet
            variable p(n,n) symmetric
            minimize   (-log_det(p))
            subject to
                p >= 0; % sdp_tol*eye(n);
                [a*p-A'*p*A  -A'*p*B 
                 -B'*p*A     (1-a)*R-B'*p*B ] >= 0; % sdp_tol*eye(2*n);
        cvx_end

        P_ellipse(:,:,i) = double(p);
        if(strcmp(cvx_status,'Solved'))
            Ellipsoid_Volume(i) = det(inv(P_ellipse(:,:,i)));
        else
            Ellipsoid_Volume(i) = NaN;
        end

    end
    [min_volume,min_index] = nanmin(Ellipsoid_Volume);
    
    % SDP Problem Data
    % p - #of outputs
    % A = A, B = -L*sqrtm(residual_var), R = eye(p)/attack_threshold
    
    A_new = A + B_original*G;
    B     = -B_original*G;
    R     = P_ellipse(:,:,min_index);
    
    % Data Structures
    a_range          = max(size(a_values));
    P_ellipse_new    = zeros(n,n,a_range);
    Ellipsoid_Volume = zeros(a_range,1);
    
    for i = 1:a_range
        clear p P_inverse
        a = a_values(i)
        
        % Solve SDP    
        cvx_begin sdp quiet
            variable p(n,n) symmetric    
            minimize   (-log_det(p))
            subject to
                p >= 0; % sdp_tol*eye(n);
                [a*p-A_new'*p*A_new  -A_new'*p*B 
                 -B'*p*A_new     (1-a)*R-B'*p*B ] >= 0; % sdp_tol*eye(2*n);             
        cvx_end
    
        P_ellipse_new(:,:,i) = double(p);    
        if(strcmp(cvx_status,'Solved'))
            Ellipsoid_Volume(i) = det(inv(P_ellipse_new(:,:,i)));
        else
            Ellipsoid_Volume(i) = NaN;
        end
    end
    
    [min_volume,min_index] = nanmin(Ellipsoid_Volume);
    PP_attack              = P_ellipse_new(:,:,min_index);
    
    % Box the output parameter
    attack_output_param.PP_attack  = PP_attack;
    attack_output_param.min_volume = min_volume;

end