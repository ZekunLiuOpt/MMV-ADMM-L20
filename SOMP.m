% Existed algorithm for solving the MMV problem based on greedy match pursuit
% Input: the measurement matrix X, the sensing matrix B, each column with unit norm, the row-sparsity K
% Output: the reconstructed sparse matrix S


function S = SOMP(X, B, K)

global BtB;
BtX = B'*X;

[numDim numBase] = size(B);  % Initialize
numSmp = size(X,2);
S = sparse(numBase, numSmp);

R = X;
IndexAtom = [];

for k = 1:K
    
    [val IndexAtom(k)] =  max(max(abs(R'*B),[],1));  % Project residual onto dictioanry and select atom as argmax
    Bk = B(:,IndexAtom);  % Orthogonal projection on the subspace spanned by all previously selected atoms
    P = pinv(Bk)*X;  % BtB(IndexAtom,IndexAtom)\BtX(IndexAtom,:);  % (Bk'*Bk)\Bk'*X;  % pinv(Bk)*X;
    % P =  invChol_mex(BtB(IndexAtom,IndexAtom))*BtX(IndexAtom,:);
    % P = inv1(BtB(IndexAtom,IndexAtom))*BtX(IndexAtom,:);
    % P = MatrixInverse(BtB(IndexAtom,IndexAtom))*BtX(IndexAtom,:);
        
    R = X - Bk*P;  % Update residual
    if norm(R,'fro')<1e-6 
        S(IndexAtom,:) = P;
        break;
    end    
end

S(IndexAtom,:) = P;
end
