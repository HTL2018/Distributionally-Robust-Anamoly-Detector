% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% MECH 6V29 - Convex Optimization in Systems & Controls
% Project - Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 8th March, 2019.
% This code is used to calculate an optimum attack threshold that will 
% result in desired false alarm rate
function output_param = compute_attack_threshold(input_param)
    
    % Problem Data    
    alarm_rate       = input_param.alarm_rate; 
    mu_attack        = input_param.mu_attack;   
    attack_cov       = input_param.attack_cov; 
    attack_tol       = 0.00001;           
    attack_high      = 1000;
    attack_low       = 0;    
    attack_threshold = 100;
    iter_counter     = 1;  
    A_matrix         = eye(size(attack_cov,1)); % attack_cov
    
    % Get Noise_Threshold by solving below SDP using bisection algorithm
    while (attack_high - attack_low > attack_tol) 
        
        iter_counter
        attack_threshold = (attack_high + attack_low)/2   
        
        cvx_begin sdp quiet
            variable Z(size(attack_cov,1),size(attack_cov,1)) symmetric
            variable z(size(attack_cov,1),1)
            variable lambda

            minimize (1 - lambda)
            subject to
                trace(A_matrix*Z) - lambda*attack_threshold >= 0;
                [Z z; z' lambda] <= [attack_cov+mu_attack*mu_attack' mu_attack; mu_attack' 1];
                [Z z; z' lambda] >= 0;
        cvx_end
        
        if(cvx_optval > 1 - alarm_rate)
            attack_high = attack_threshold;            
        elseif(cvx_optval < 1 - alarm_rate)
            attack_low  = attack_threshold;             
        end
               
        threshold(iter_counter) = attack_threshold;
        iter_counter            = iter_counter + 1;

    end            
    
    output_param.attack_threshold = attack_threshold;  
    output_param.z                = z;
    output_param.Z                = Z;
    output_param.lambda           = lambda;
    
end