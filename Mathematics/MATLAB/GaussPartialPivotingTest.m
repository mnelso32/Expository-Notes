function [x] = GaussPartialPivotingTest(A,b)

n = length(b);

for k=1:n-1    
    % First, find the row to pivot with    

    % option a: (use max())    
    % [pivotVal, pivotRow] = max(abs(A(k:n,k)));    
    % adjust to row number in A (there are k-1 above the diagonal)   
    % pivotRow = pivotRow+k-1;
    
    % option b: (manually)    
    pivotRow=k;    
    for j=k+1:n       
        if abs(A(j,k)) > abs(A(pivotRow,k))           
            pivotRow=j;        
        end    
    end
    
    % Interchange rows of matrix and RHS vector if pivot is not row k
    if (pivotRow ~= k)
        A([k pivotRow],:) = A([pivotRow k],:);
        b([k pivotRow]) = b([pivotRow k]);
    end
    
    % Apply transformation to remaining submatrix and RHS vector
    % (same as in GaussElimination)
    for i=k+1:n
        m = A(i,k)/A(k,k); % multipliyer for current row i
        for j=k:n % loop to update row i
            A(i,j) = A(i,j) - m*A(k,j);
        end
        b(i) = b(i) - m*b(k); % update RHS
    end
end

A
b