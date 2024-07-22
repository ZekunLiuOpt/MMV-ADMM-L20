% Function that solves the L20-projection subproblem
% Input: a given matrix C, with the known number of non-zero rows s
% Output: the L20-projection matrix T
% Written by: Zekun Liu (11/02/2023)
% Reference:
% [1] Z. Liu and S. Yu. Alternating Direction Method of Multipliers Based on $\ell_{2,0}$-Norm for Multiple Measurement Vector Problem.
%     IEEE Transactions on Signal Processing, vol. 71, pp. 3490-3501, 2023.
% Latest Revision: 22/07/2024


function T = rowshrinkL20(C, s)

m = size(C, 1);
n = size(C, 2);
T = zeros(m, n);

L2 = [];
for i = 1:m
    L2(i) = norm(C(i, :), 2);
end

[~, index] = sort(L2, 'descend');

for i = 1:s
    T(index(i), :) = C(index(i), :);
end
end
