% Our algorithm for solving the MMV problem
% Input: the measurement matrix Y, the sensing matrix Phi, the row-sparsity s (if one has the prior information; otherwise using the formula in the paper to calculate s first, and one can notice that 
% spark(Phi)=rank(Phi)+1=M+1 holds with high probability), the penalty parameter rho>0
% Output: the reconstructed sparse matrix S
% Copyright: Zekun Liu


function S = MMV_ADMM_L20(Y,Phi,s,rho)

c = size(Y,2);  % Get dimensions of S
r = size(Phi,2); 

L = zeros(r,c);  % Initialize Lagragian to be nothing (seems to work well)
maxIter = 1000;  % Set the maximum number of iterations (make really big to ensure convergence)
I = speye(r);  % Set the sparse identity matrix
S = randn(r,c);  % Initialize S randomly

epsilon1 = 1e-6;  %stopping tolerance
epsilon2 = 1e-6;
epsilon3 = 1e-6;

invA = inv(2*Phi'*Phi+rho*I);  % calculate the inverse only once outside the iterations

for n = 1:maxIter
    Sold = S;
    C = rowshrinkL20(S - (L/rho),s);
    S = invA*(2*Phi'*Y + rho*C + L); 
    L = L + rho*(C - S);  

    rp = norm(S-C,'fro');  % Check if satisfies the expected tolerance. For comparisons, one can also only check if rd<epsilon2
    rd = norm(Sold-S,'fro');
    rl = norm(L,'fro');
    if  (rp < epsilon1)  && (rd < epsilon2)  &&  (rl < epsilon3) 
        break;
    end
end

S = C;
end
