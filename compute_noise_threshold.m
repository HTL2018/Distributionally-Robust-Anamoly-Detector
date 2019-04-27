% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 13th March, 2019.
% This code is used to calculate an optimum noise threshold that will 
% result in desired false alarm rate

function out_param = compute_noise_threshold(input_param)

    % Problem Data
    
    alarm_rate      = input_param.alarm_rate; 
    mu_noise        = input_param.mu_noise;   
    Sigma_noise     = input_param.Sigma_noise; 
    noise_tol       = 0.00001;           
    noise_high      = 100;
    noise_low       = 0;    
    noise_threshold = 100;
    iter_counter    = 1;    
    
    % Get Noise_Threshold by solving below SDP using bisection algorithm
    while (noise_high - noise_low > noise_tol)         
        
        noise_threshold = (noise_high + noise_low)/2
        
        cvx_begin sdp quiet
            variable Z(size(Sigma_noise,1),size(Sigma_noise,1)) symmetric
            variable z(size(Sigma_noise,1),1)
            variable lambda

            minimize (1 - lambda)
            subject to
                trace(inv(Sigma_noise)*Z) - lambda*noise_threshold >= 0;
                [Z z; z' lambda] <= [Sigma_noise+mu_noise*mu_noise' mu_noise; mu_noise' 1];
                [Z z; z' lambda] >= 0;
        cvx_end
        
        if(cvx_optval > 1 - alarm_rate)
            noise_high = noise_threshold;            
        elseif(cvx_optval < 1 - alarm_rate)
            noise_low  = noise_threshold;              
        end
        
        threshold(iter_counter) = noise_threshold;
        iter_counter            = iter_counter + 1;

    end        
    
    out_param.noise_threshold = noise_threshold;  
    out_param.z               = z;
    out_param.Z               = Z;
    out_param.lambda          = lambda;

end