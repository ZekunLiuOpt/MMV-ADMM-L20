% Function that solves the L20-projection subproblem
% Input: a given matrix C, with the known number of non-zero rows s
% Output: the L20-projection matrix T
% Written by: Zekun Liu


function T = rowshrinkL20(C,s)

m = size(C,1);
n = size(C,2);
T = zeros(m,n);

L2 = [];
for i = 1:m
    L2(i) = norm(C(i,:),2);
end

[~,index] = sort(L2,'descend');

for i = 1:s
    T(index(i),:) = C(index(i),:);
end
end
