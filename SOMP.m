% Existed algorithm for solving the MMV problem based on greedy match pursuit
%
% Input: 
%         X: the measurement matrix
%         B: the sensing matrix, each column with unit norm
%         K: the row-sparsity K
%
% Output: 
%         S: the reconstructed sparse matrix
%
% Written by jk123vip (https://github.com/jk123vip/SRC_SOMP_matlab/blob/master/SOMP.m)
%
% Reference: 
% [1] S. F. Cotter, B. D. Rao, K. Engan, and K. Kreutz-Delgado. 
%     Sparse solutions to linear inverse problems with multiple measurement vectors.
%     IEEE Trans. Signal Process., vol. 53, no. 7, pp. 2477–2488, 2005.
% [2] J. A. Tropp, A. C. Gilbert, and M. Strauss. 
%     Algorithms for simultaneous sparse approximation. Part I: Greedy pursuit.
%     Signal Process., vol. 86, no. 3, pp. 572–588, Mar. 2006.
%
% Latest Revision: 17/10/2024


function S = SOMP(X, B, K)

global BtB;
BtX = B' * X;

% initialization
[numDim numBase] = size(B);  
numSmp = size(X, 2);
S = sparse(numBase, numSmp);

R = X;
IndexAtom = [];

for k = 1:K
    
    [val IndexAtom(k)] =  max(max(abs(R' * B), [], 1));  % project residual onto dictioanry and select atom as argmax
    Bk = B(:, IndexAtom);  % orthogonal projection on the subspace spanned by all previously selected atoms
    P = pinv(Bk) * X;  % BtB(IndexAtom, IndexAtom) \ BtX(IndexAtom, :);  % (Bk' * Bk) \ Bk' * X;  % pinv(Bk) * X;
    % P = invChol_mex(BtB(IndexAtom, IndexAtom)) * BtX(IndexAtom, :);
    % P = inv1(BtB(IndexAtom, IndexAtom)) * BtX(IndexAtom, :);
    % P = MatrixInverse(BtB(IndexAtom, IndexAtom)) * BtX(IndexAtom, :);
        
    R = X - Bk * P;  % update residual
    if norm(R, 'fro') < 1e-6 
        S(IndexAtom, :) = P;
        break;
    end    
end

S(IndexAtom, :) = P;
end

S(IndexAtom, :) = P;
end
