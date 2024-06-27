% Function that solves the L21-shrink subproblem argmin_{T} kaka*||T||_{2,1}+1/2*||T-C||^2
% Input: a given matrix C, with the parameter kaka
% Output: the L21-shrink solution matrix T
% Written by: Zekun Liu (10/02/2023)


function T = rowshrinkL21(C,kaka)

m = size(C,1);
n = size(C,2);
T = zeros(m,n);

for i = 1:m
    l = norm(C(i,:),2);
    if l > kaka
        T(i,:) = (l-kaka)/l*C(i,:);
    end
end
end
