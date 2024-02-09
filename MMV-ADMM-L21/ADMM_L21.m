% Existed algorithm for solving the MMV problem based on L21-norm
% Input: the measurement matrix X, the sensing matrix A, the parameter gamma>0, rho>0
% Output: the reconstructed sparse matrix B
% Copyright: Zekun Liu


function B = ADMML21(X, A, gamma, rho)

c = size(X,2);  % Get dimensions of B
r = size(A,2); 

L = zeros(r,c);  % Initialize Lagragian to be nothing (seems to work well)
maxIter = 1000;  % Set the maximum number of iterations (make really big to ensure convergence)
I = speye(r);  % Set the sparse identity matrix
C = randn(r,c);  % Initialize C randomly
epsilon1 = 1e-6;  %stopping criteria 
epsilon2 = 1e-6;

invA = inv(A'*A+rho*I);

for n = 1:maxIter
    Cold = C;
    B = invA*(A'*X + rho*C - L);
    C = rowshrinkL21(B + (L/rho),gamma/rho); 
    L = L + rho*(B - C);  

    rp = norm(B-C,'fro');
    rd = norm(Cold-C,'fro');
    if (rp < epsilon1) && (rd < epsilon2)
        break;
    end
end

B = C;
end
