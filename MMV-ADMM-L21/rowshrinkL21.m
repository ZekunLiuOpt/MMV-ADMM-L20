% Function that solves the L21-shrink subproblem argmin_{T} kaka*||T||_{2,1}+1/2*||T-C||^2
% Input: a given matrix C, with the parameter kaka
% Output: the L21-shrink solution matrix T
% Written by: Zekun Liu (10/02/2023)
% Reference:
% [1] Q. Qu, N. M. Nasrabadi, and T. D. Tran. Abundance estimation for bilinear mixture models via joint sparse 
%     and low-rank representation. IEEE Trans. Geosci. Remote Sens., vol. 52, no. 7, pp. 4404–4423, Jul.2014.
% Latest Revision: 23/07/2024


function T = rowshrinkL21(C, kaka)

m = size(C, 1);
n = size(C, 2);
T = zeros(m, n);

for i = 1:m
    l = norm(C(i, :), 2);
    if l > kaka
        T(i, :) = (l - kaka) / l * C(i, :);
    end
end
end
