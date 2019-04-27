% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu 
% Date: 1st March, 2019.

function attack_output_param = attack_bounding_ellipsoid(attack_input_param)

    % problem Data
    A            = attack_input_param.A;
    B            = attack_input_param.B;
    L            = attack_input_param.L;
    K            = attack_input_param.K;
    B_original   = B;
    mu_attack    = attack_input_param.mu_attack;
    alarm_rate   = attack_input_param.alarm_rate;
    attack_cov   = attack_input_param.attack_cov;
    residual_var = attack_input_param.residual_var;
    a_values     = 0.01:0.1:0.99;    
    n            = size(A,1); 
    type         = attack_input_param.type;
    
    if type == 1 
        % Type = 1 - DR Case
        
        attack_threshold = attack_input_param.attack_threshold;
    else
        % Type = 2 - Chi Squared Case
        attack_threshold = ncx2inv(1-alarm_rate,n,0);
    end   
    
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
    
    A_new = A + B_original*K;
    B     = -B_original*K;
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
    % Box the output parameter
    attack_output_param.P          = P_ellipse_new(:,:,min_index);
    attack_output_param.min_volume = min_volume;

end