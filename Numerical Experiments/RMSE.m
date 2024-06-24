% Function that calculates the RMSE 
% Input: a recovered matrix x1, a ground truth x2
% Output: the RMSE of x1 to x2
% Written by: Zekun Liu


function y = RMSE(x1,x2)

N = size(x1,1);
J = size(x1,2);

y = norm(x1-x2,'fro')/sqrt(N*J);
end
