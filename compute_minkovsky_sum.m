% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 13th March, 2019.

function minkov_sum_matrix = compute_minkovsky_sum(matrix_1,matrix_2)
% Given two matrices correspondin to two ellipsoids, this function computes
% the minkovsky sum of both ellipsoids and returns the inverse of the matrix
% associated with the minkovsky sum.

    inv_matrix_1      = inv(matrix_1); 
    inv_matrix_2      = inv(matrix_2);    
    minkov_sum_matrix = inv((sqrt(trace(inv_matrix_1)) + sqrt(trace(inv_matrix_2)))*...
                        ((inv_matrix_1/sqrt(trace(inv_matrix_1))) + ...
                        (inv_matrix_2/ sqrt(trace(inv_matrix_2)))));

end

