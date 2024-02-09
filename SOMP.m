function S = SOMP(X, B, K)
% 多源数据DCS的同步正交匹配追踪算法
% Input:
% X - samples matrix, one sample each column, with unit norm
% B - dictionary, each column with unit norm
% K - spasity degree
% Output:
% S - sparse coding matrix

% 

global BtB;

BtX = B'*X;% Ô¤¼ÆËã
% initial

[numDim numBase] = size(B);
numSmp = size(X,2);
S = sparse(numBase, numSmp);

R = X;
IndexAtom = [];

% main loop

for k = 1:K
    
    % project residual onto dictioanry and select atom as arg max
    [val IndexAtom(k)] =  max(max(abs(R'*B),[],1)); % inf norm

    % Orthogonal projection on the subspace spanned by all previously selected atoms

    Bk = B(:,IndexAtom);
    P = pinv(Bk)*X;%BtB(IndexAtom,IndexAtom)\BtX(IndexAtom,:); % (Bk'*Bk)\Bk'*X; % pinv(Bk)*X;%
%P =  invChol_mex(BtB(IndexAtom,IndexAtom))*BtX(IndexAtom,:);
% P = inv1(BtB(IndexAtom,IndexAtom))*BtX(IndexAtom,:);
% P = MatrixInverse(BtB(IndexAtom,IndexAtom))*BtX(IndexAtom,:);
    
    % update residual
    
    R = X - Bk*P;
    
    if norm(R,'fro')<1e-6%总残差足够小直接结束,注意此部分仅在K不确定是否是实际稀疏支撑集稀疏度时使用
        %如果确定一致此部分不需要加，否则增加运行时间
        S(IndexAtom,:) = P;
        break;
    end
    
end

S(IndexAtom,:) = P;

end
