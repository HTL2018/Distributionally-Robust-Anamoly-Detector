% Venkatraman Renganathan, Navid Hashemi
% Email: (vrengana, navid.hashemi)@utdallas.edu
% Distributionally Robust Ellipsoidal Bounds for Reachable Sets
% Date: 13th March, 2019.

function minkov_sum_matrix = compute_minkovsky_sum_Ver1(matrix_1,matrix_2)
% Given two matrices correspondin to two ellipsoids, this function computes
% the minkovsky sum of both ellipsoids and returns the matrix associated
% with the minkovsky sum.
    
    minkov_sum_matrix = (sqrt(trace(matrix_1)) + sqrt(trace(matrix_2)))*...
                        ((matrix_1/sqrt(trace(matrix_1))) + ...
                        (matrix_2/ sqrt(trace(matrix_2))));

end

