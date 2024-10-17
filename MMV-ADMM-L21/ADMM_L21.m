% Existed algorithm for solving the MMV problem based on L21-norm
%
% Input: 
%         X: the measurement matrix
%         A: the sensing matrix
%     gamma: the regularization paramter (> 0)
%       rho: the penalty parameter (> 0)
%
% Output: 
%         B: the reconstructed sparse matrix
%
% Written by Zekun Liu, 10/02/2023
%
% Reference:
% [1] Q. Qu, N. M. Nasrabadi, and T. D. Tran. 
%     Abundance estimation for bilinear mixture models via joint sparse and low-rank representation.
%     IEEE Trans. Geosci. Remote Sens., vol. 52, no. 7, pp. 4404–4423, Jul.2014.
%
% Latest Revision: 17/10/2024


function B = ADMM_L21(X, A, gamma, rho)

c = size(X, 2);   % Get the dimension of B
r = size(A, 2); 

L = zeros(r, c);  % Initialize Lagragian to be nothing (seems to work well)
maxIter = 1000;   % Set the maximum number of iterations (make really big to ensure convergence)
I = speye(r);     % Set the sparse identity matrix
C = randn(r, c);  % Initialize C randomly
epsilon1 = 1e-6;  % Stopping criteria 
epsilon2 = 1e-6;

invA = inv(A' * A + rho * I);

for n = 1:maxIter
    Cold = C;
    B = invA * (A' * X + rho * C - L);
    C = rowshrinkL21(B + (L / rho), gamma / rho); 
    L = L + rho * (B - C);  

    % For comparisons, one can also only check if rd < epsilon2
    rp = norm(B - C, 'fro');
    rd = norm(Cold - C, 'fro');
    if (rp < epsilon1) && (rd < epsilon2)  
    
        break;
    end
end

B = C;
end
