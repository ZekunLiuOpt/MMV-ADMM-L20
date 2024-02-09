function Xnew = sniht(Y,A,k,X0supp,printitn)
% 
% Simultaneous Normalized Iterative Hard Thresholding algorihtm proposed 
% in Blanchard et al (2014). 
% Single measurement vector case is NIHT in Blumensath and Davies (2010). 
% 
%  INPUT  
%        Y  :=  matrix of measurements (response)
%        A  :=  matrix of predictors (measurement matrix) 
%        k  :=  the number of nonzero coefficients 
%        printitn := print iteration number (modulo)   
% OUTPUT  
%        Xnew    := estimated signal parameter matrix with k nonzero rows
%        X1supp  := estimated support set of non-zeros
%        failure := equal to 1 if algorithm failed
%
% Author: E. Ollila, Oct 15th, 2014
%---------------------------------------------------


[m, n] = size(A);
[~, q] = size(Y);

if nargin < 5, printitn=0; end
if nargin < 4, X0supp=[]; end


%%-- initial approximation is the zero matrix 
X0  = zeros(n,q); 

ARRAYSIM =  false;

%%- DetectSupport
if isempty(X0supp),
    R = A'*Y;
    [~, indx] = sort(sum(R.*conj(R),2),'descend');
    X0supp = indx(1:k);
    
    if ARRAYSIM,
        if abs(indx(2) - indx(1)) ==1, 
             X0supp(2)=indx(3); 
        end 
     
        if  abs(X0supp(2) - indx(1)) ==1,
            X0supp(2)=indx(4); 
        end
        %X0supp = sort(X0supp);
        %X0supp'
    end
end

ITERMAX = 1000;
%cc = 0.02;
objold = 1e12; % a big number
failure = 0;

for iter = 1:ITERMAX
    
    %%-- Compute the negative gradient 
    R = Y - A*X0;
    G = A'*R; 

    %%- stepsize
    mu = norm(G(X0supp,:),'fro')/norm(A(:,X0supp)*G(X0supp,:),'fro');
    mu = mu^2;

    %%-- Next proposal for b 
    X1 = X0 + mu*G;
    %%-- Detect Support
    [~, indx] = sort(sum(X1.*conj(X1),2),'descend');
    X1supp = indx(1:k);
    
    %%-- Threshold
    Xnew   = zeros(n,q);
    Xnew(X1supp,:)= X1(X1supp,:);

    %%-- Stopping criteria          
    crit = norm(Xnew-X0,'fro')^2/norm(Xnew,'fro')^2;
    %crit = norm(bnew-b0)^2;
       
    objnew =   sum(sum((Y-A*Xnew).*conj(Y-A*Xnew)));  % = norm(Y-A*Xnew,'fro')^2
    
    if mod(iter,printitn) == 0
       fprintf('iter %2d: mu = %f objnew = %f crit = %f\n',iter,mu,objnew,crit);  
    end
    
    if objnew > objold 
        failure = 1;
    end
    
    if crit < 1e-14
      break 
    end 
        
    X0   = Xnew;
    X0supp = X1supp;
        
end

%fprintf('NORIHT terminating at iter = %d crit = %f\n',iter,crit)  

if iter == ITERMAX, failure=1; end;

X1supp = sort(X1supp);
