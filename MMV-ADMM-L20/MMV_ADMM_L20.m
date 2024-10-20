% Our algorithm for solving the MMV problem based on L20-norm
%
% Input: 
%         Y: the measurement matrix
%       Phi: the sensing matrix
%         s: the row-sparsity
%            directly set it if one has the prior information
%            otherwise using the formula in the paper to calculate it first
%            one can notice that spark(Phi) = rank(Phi) + 1 = M + 1 holds 
%            with high probability for Gaussian random matrix
%       rho: the penalty parameter (> 0)
%
% Output: 
%         S: the reconstructed sparse matrix
%
% Written by Zekun Liu, 20/02/2023
%
% Reference:
% [1] Z. Liu and S. Yu. 
%     Alternating Direction Method of Multipliers Based on $\ell_{2,0}$-Norm for Multiple Measurement Vector Problem.
%     IEEE Transactions on Signal Processing, vol. 71, pp. 3490-3501, 2023.
%
% Latest Revision: 17/10/2024


function S = MMV_ADMM_L20(Y, Phi, s, rho)

c = size(Y, 2);   % Get the dimension of S
r = size(Phi, 2); 

L = zeros(r, c);  % Initialize Lagragian to be nothing (seems to work well)
maxIter = 1000;   % Set the maximum number of iterations (make really big to ensure convergence)
I = speye(r);     % Set the sparse identity matrix
S = randn(r, c);  % Initialize S randomly

epsilon1 = 1e-6;  % Stopping tolerance
epsilon2 = 1e-6;
epsilon3 = 1e-6;

invA = inv(2 * Phi' * Phi + rho * I);  % calculate the inverse only once outside the iterations
% Use the SMW-formula to calculate the inverse, and call it MMV-ADMM-L20-SMW 
% invA = I / rho - 2 * Phi' / (eye(M) + 2 * Phi * Phi' / rho) * Phi / (rho^2); 

for n = 1:maxIter
    Sold = S;
    C = rowshrinkL20(S - (L / rho), s);
    S = invA * (2 * Phi' * Y + rho * C + L); 
    L = L + rho * (C - S);  

    % Check if stop criterions satisfy the expected tolerance 
    % For comparisons, one can also only check if rd < epsilon2
    % Call it MMV-ADMM-L20-NCC when skips the stop check
    % Call it MMV-ADMM-L20-SeeCC when iterates to MaxIter and records all of the rp, rd, and rl
    rp = norm(S - C, 'fro'); 
    rd = norm(Sold - S, 'fro');
    rl = norm(L, 'fro');
    if  (rp < epsilon1)  && (rd < epsilon2)  &&  (rl < epsilon3) 
        break;
    end
end

S = C;
end
