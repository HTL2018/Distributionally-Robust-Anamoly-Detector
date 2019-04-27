% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 15th March, 2019.
% This code is used to calculate an optimum residual threshold that will 
% result in desired false alarm rate

function residual_threshold = compute_residual_threshold(input_param)

    % Problem Data
    
    alarm_rate         = input_param.alarm_rate; 
    mu_residual        = input_param.mu_residual;   
    Sigma_residual     = input_param.Sigma_residual; 
    residual_tol       = 0.00001;           
    residual_high      = 100;
    residual_low       = 0;    
    residual_threshold = 100;
    iter_counter       = 1;    
    
    % Get Noise_Threshold by solving below SDP using bisection algorithm
    while (residual_high - residual_low > residual_tol)         
        
        residual_threshold = (residual_high + residual_low)/2
        
        cvx_begin sdp quiet
            variable Z(size(Sigma_residual,1),size(Sigma_residual,1)) symmetric
            variable z(size(Sigma_residual,1),1)
            variable lambda

            minimize (1 - lambda)
            subject to
                trace(inv(Sigma_residual)*Z) - lambda*residual_threshold >= 0;
                [Z z; z' lambda] <= [Sigma_residual+mu_residual*mu_residual' mu_residual; mu_residual' 1];
                [Z z; z' lambda] >= 0;
        cvx_end
        
        if(cvx_optval > 1 - alarm_rate)
            residual_high = residual_threshold;            
        elseif(cvx_optval < 1 - alarm_rate)
            residual_low  = residual_threshold;              
        end
        
        threshold(iter_counter) = residual_threshold;
        iter_counter            = iter_counter + 1;

    end        

end