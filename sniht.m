% Existed algorithm for solving the MMV problem based on IHT
%
% Input: 
%         Y: the measurement matrix Y
%         A: the sensing matrix
%         k: the row-sparsity
%
% Output: 
%      Xnew: the reconstructed sparse matrix
%
% Written by E. Ollila, 15/10/2014 (https://github.com/AmmarMian/huber_mm_framework/
% blob/9b2f84772e21807c6ea6fb72c18d33c99b51b2e6/python/notebooks/matlab/sniht.m#L4)
%
% Reference: 
% [1] J. D. Blanchard, M. Cermak, D. Hanle, and Y. Jing.
%     Greedy algorithms for joint sparse recovery. 
%     IEEE Trans. Signal Process., vol. 62, no. 7, pp. 1694–1704, Apr. 2014.
%
% Latest Revision: 17/10/2024


function Xnew = sniht(Y, A, k, X0supp, printitn) 

[m, n] = size(A);
[~, q] = size(Y);

if nargin < 5, printitn = 0; end
if nargin < 4, X0supp = []; end

% initial approximation is the zero matrix 
X0  = zeros(n, q); 

ARRAYSIM =  false;

% DetectSupport
if isempty(X0supp)
    R = A' * Y;
    [~, indx] = sort(sum(R .* conj(R), 2), 'descend');
    X0supp = indx(1:k);
    
    if ARRAYSIM
        if abs(indx(2) - indx(1)) == 1 
             X0supp(2) = indx(3); 
        end 
     
        if  abs(X0supp(2) - indx(1)) == 1
            X0supp(2) = indx(4); 
        end
        % X0supp = sort(X0supp);
        % X0supp'
    end
end

ITERMAX = 1000;
% cc = 0.02;
objold = 1e12;  % a big number
failure = 0;

for iter = 1:ITERMAX
    
    % Compute the negative gradient 
    R = Y - A * X0;
    G = A' * R; 

    % stepsize
    mu = norm(G(X0supp, :), 'fro') / norm(A(:, X0supp) * G(X0supp, :), 'fro');
    mu = mu^2;

    % Next proposal for b 
    X1 = X0 + mu * G;
    % Detect Support
    [~, indx] = sort(sum(X1 .* conj(X1), 2), 'descend');
    X1supp = indx(1:k);
    
    % Threshold
    Xnew = zeros(n, q);
    Xnew(X1supp, :) = X1(X1supp, :);

    % Stopping criteria          
    crit = norm(Xnew - X0, 'fro')^2 / norm(Xnew, 'fro')^2;
    % crit = norm(bnew - b0)^2;
       
    objnew = sum(sum((Y - A * Xnew) .* conj(Y - A * Xnew)));  % = norm(Y - A * Xnew, 'fro')^2
    
    if mod(iter, printitn) == 0
       fprintf('iter %2d: mu = %f objnew = %f crit = %f\n', iter, mu, objnew, crit);  
    end
    
    if objnew > objold 
        failure = 1;
    end
    
    if crit < 1e-14
      break 
    end 
        
    X0 = Xnew;
    X0supp = X1supp;
        
end

% fprintf('NORIHT terminating at iter = %d crit = %f\n', iter, crit)  
if iter == ITERMAX, failure = 1; end
X1supp = sort(X1supp);
end
