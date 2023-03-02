function [x] = BandedGaussElimination(A,b,k)

n = length(b);

for j=1:n-1     
    % Check to see if the pivot is zero    
    if abs(A(j,j)) < 1e-15        
        error('A has diagonal entries of zero')    
    end    
    
    % Apply transformation to remaining submatrix and RHS vector
    for i=j+1:j+k
        m = A(i,j)/A(j,j); % multipliyer for current row i
        for j1=j:j+k % loop to update row i
            A(i,j1) = A(i,j1) - m*A(j,j1);
        end
        b(i) = b(i) - m*b(j); % update RHS
    end
end

% A is now upper triangular. Use backsubstitution to solve the transformed problem.
x = BackSubstitution(A,b);