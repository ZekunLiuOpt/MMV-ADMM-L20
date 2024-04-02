% Experiments that comparing the improvement with SMW-formula to MMV-ADMM-L20
% Copyright: Zekun Liu


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

    Index_K = randperm(N);
    spar_arr = randn(1,K*J);
    S = sparse(repelem(Index_K(1:K),J),repmat([1:1:J],1,K),spar_arr,N,J);
    Phi = sqrt(1/M)*randn(M,N);
    Y = Phi*S;
    
    tic
    S1 = MMV_ADMM_L20(Y, Phi, K, 1);
    toc

    timesaver1 = [timesaver1,toc];
    RMSEsaver1 = [RMSEsaver1,RMSE(S1,S)];

    tic
    S2 = MMV_ADMM_L20_SMW(Y, Phi, K, 1);
    toc

    timesaver2 = [timesaver2,toc];
    RMSEsaver2 = [RMSEsaver2,RMSE(S2,S)];
end

mean(RMSEsaver1)
std(RMSEsaver1)
mean(RMSEsaver2)
std(RMSEsaver2)

mean(timesaver1)
std(timesaver1)
mean(timesaver2)
std(timesaver2)

save test7
load test7
