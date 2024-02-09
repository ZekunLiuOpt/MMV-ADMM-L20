clear all
clc

%Set up
N = 5000;
M = 1500;
K = 500;
J = 10;

timesaver1 = [];
RMSEsaver1 = [];
timesaver2 = [];
RMSEsaver2 = [];

T = 10;

for i = 1:T

    Index_K = randperm(N);%Sparse support set
    
    spar_arr = randn(1,K*J);%The value of sparse matrix to be recovered
    
    S = sparse(repelem(Index_K(1:K),J),repmat([1:1:J],1,K),spar_arr,N,J);%N*J sparse matrix to be recovered
    
    Phi = sqrt(1/M)*randn(M,N);%M*N unit-norm Gaussian random matrix
    
    Y = Phi*S;%Measurement value matrix
    
    tic
    S1 = MMV_ADMM_L20(Y, Phi, K, 1);%S1 is the recovered estimation of S.
    toc

    timesaver1 = [timesaver1,toc];

    RMSEsaver1 = [RMSEsaver1,RMSE(S1,S)];

    tic
    S2 = MMV_ADMM_L20_SMW(Y, Phi, K, 1);%S1 is the recovered estimation of S.
    toc

    timesaver2 = [timesaver2,toc];

    RMSEsaver2 = [RMSEsaver2,RMSE(S2,S)];

end

save test9
load test9

mean(RMSEsaver1)
std(RMSEsaver1)
mean(RMSEsaver2)
std(RMSEsaver2)

mean(timesaver1)
std(timesaver1)
mean(timesaver2)
std(timesaver2)
