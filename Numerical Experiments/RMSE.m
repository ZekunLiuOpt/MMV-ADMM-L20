% Function that calculates the RMSE 
%
% Input: 
%         x1: the recovered matrix
%         x2: the ground truth
%
% Output: 
%          y: the RMSE of x1 and x2
%
% Written by Zekun Liu, 02/03/2023
%
% Latest Revision: 17/10/2024


function y = RMSE(x1, x2)

N = size(x1, 1);
J = size(x1, 2);

y = norm(x1 - x2, 'fro') / sqrt(N * J);
end
