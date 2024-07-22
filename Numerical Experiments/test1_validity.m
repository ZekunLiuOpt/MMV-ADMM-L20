% Experiments that test the validity of our algorithm, applying MMV-ADMM-L20-NCC (without stop criterion)
% Written by: Zekun Liu (03/03/2023)
% Latest Revision: 22/07/2024


clear all
clc

% Set up
N = 500;
M = 150;
K = 50;
J = 10;

timesaver = [];
RMSEsaver = [];

T = 10;
for i = 1:T

    Index_K = randperm(N);             % Sparse support set
    spar_arr = randn(1, K * J);        % The value of sparse matrix to be recovered
    S = sparse(repelem(Index_K(1:K), J), repmat([1:1:J], 1, K), spar_arr, N, J);   % N × J sparse matrix to be recovered
    Phi = sqrt(1 / M) * randn(M, N);   % M × N unit-norm Gaussian random matrix
    Y = Phi * S;                       % Measurement matrix
    
    tic
    S1 = MMV_ADMM_L20_NCC(Y, Phi, K, 1);  % S1 is the recovered estimation of S.
    toc
    
    timesaver = [timesaver, toc];
    RMSEsaver = [RMSEsaver, RMSE(S1, S)];
end

mean(RMSEsaver)
std(RMSEsaver)
mean(timesaver)
std(timesaver)

save test1
load test1
