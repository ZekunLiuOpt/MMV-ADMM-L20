% Experiments that test the convergence criterion changes averaged over T trials, applying MMV-ADMM-L20-SeeCC (record stop criterion)
% Copyright: Zekun Liu


clear all
clc

%Set up
N = 500;
M = 150;
K = 50;
J = 10;

criterion1saver=[];
criterion2saver=[];
criterion3saver=[];

T = 10;
for i = 1:T

    Index_K = randperm(N);
    spar_arr = randn(1,K*J);    
    S = sparse(repelem(Index_K(1:K),J),repmat([1:1:J],1,K),spar_arr,N,J);
    Phi = sqrt(1/M)*randn(M,N);
    Y = Phi*S;
    
    tic
    [S1,c1,c2,c3] = MMV_ADMM_L20_SeeCC(Y, Phi, K, 1);
    toc

    criterion1saver = [criterion1saver;c1];
    criterion2saver = [criterion2saver;c2];
    criterion3saver = [criterion3saver;c3];
end

Ave_C1 = mean(criterion1saver);
Ave_C2 = mean(criterion2saver);
Ave_C3 = mean(criterion3saver);

%plot the average change of convergence criterion
figure;
ind = 1:1000;
plot(ind,Ave_C1,'LineWidth',1.25);
xlim([0 1000]);
xlabel('iteration');
ylabel('norm(B^{k}-S^{k},''fro'')');
title('The change of the first convergence criterion in iteration');

figure;
ind = 1:1000;
plot(ind,Ave_C2,'LineWidth',1.25);
xlim([0 1000]);
xlabel('iteration');
ylabel('norm(S^{k+1}-S^{k},''fro'')');
title('The change of the second convergence criterion in iteration');

figure;
ind = 1:1000;
plot(ind,Ave_C3,'LineWidth',1.25);
xlim([0 1000]);
xlabel('iteration');
ylabel('norm(L^{k},''fro'')');
title('The change of the third convergence criterion in iteration');

%compare the influence of the convergence criterion

timesaver1 = [];  %No convergence criterion
RMSEsaver1 = [];
timesaver2 = [];  %With convergence criterion
RMSEsaver2 = [];

for i = 1:T
    Index_K = randperm(N);
    spar_arr = randn(1,K*J);
    S = sparse(repelem(Index_K(1:K),J),repmat([1:1:J],1,K),spar_arr,N,J);
    Phi = sqrt(1/M)*randn(M,N);
    Y = Phi*S;
    
    tic
    S1 = MMV_ADMM_L20_NCC(Y, Phi, K, 1);
    toc

    timesaver1 = [timesaver1,toc];
    RMSEsaver1 = [RMSEsaver1,RMSE(S1,S)];

    tic
    S2 = MMV_ADMM_L20(Y, Phi, K, 1);
    toc

    timesaver2 = [timesaver2,toc];
    RMSEsaver2 = [RMSEsaver2,RMSE(S2,S)];
end

mean(RMSEsaver1)
mean(RMSEsaver2)
std(RMSEsaver1)
std(RMSEsaver2)

mean(timesaver1)
mean(timesaver2)
std(timesaver1)
std(timesaver2)

%plot the average change of convergence criterion in 200 iterations
figure;
ind = 1:200;
plot(ind,Ave_C1(1:200),'LineWidth',1.25);
xlim([0 200]);
xlabel('iteration','FontSize',15);
ylabel('$\left \| B^{k}-S^{k} \right \|_{F}$','interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',17);
%title('The average change of the first convergence criterion','FontSize',15);

figure;
ind = 1:200;
plot(ind,Ave_C2(1:200),'LineWidth',1.25);
xlim([0 200]);
xlabel('iteration','FontSize',15);
ylabel('$\left \| S^{k+1}-S^{k} \right \|_{F}$','interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',17);
%title('The average change of the second convergence criterion','FontSize',15);

figure;
ind = 1:200;
plot(ind,Ave_C3(1:200),'LineWidth',1.25);
xlim([0 200]);
xlabel('iteration','FontSize',15);
ylabel('$\left \| L^{k} \right \|_{F}$','interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',17);
%title('The average change of the third convergence criterion','FontSize',15);

save test2
load test2
