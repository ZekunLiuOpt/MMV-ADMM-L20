% Experiments that comparing the performance of different algorithms under different sensor numbers
% Copyright: Zekun Liu


clear all
clc

%Set up
N = 500;
M = 150;
K = 50;
options = spgSetParms('iterations',1000,'bpTol',1e-6,'optTol',1e-6);

J_set = [1:20];

Timesaver1 = [];
Timesaver2 = [];
Timesaver3 = [];
Timesaver4 = [];
Timesaver5 = [];
Timesaver6 = [];

RMSEsaver1 = [];
RMSEsaver2 = [];
RMSEsaver3 = [];
RMSEsaver4 = [];
RMSEsaver5 = [];
RMSEsaver6 = [];

CNT = 100;
for i = 1:length(J_set)
    J = J_set(i);

    timesaver1 = [];
    timesaver2 = [];
    timesaver3 = [];
    timesaver4 = [];
    timesaver5 = [];
    timesaver6 = [];
    
    rmsesaver1 = [];
    rmsesaver2 = [];
    rmsesaver3 = [];
    rmsesaver4 = [];
    rmsesaver5 = [];
    rmsesaver6 = [];
    
    for j = 1:CNT
        Index_K = randperm(N);
        spar_arr = randn(1,K*J);
        S = sparse(repelem(Index_K(1:K),J),repmat([1:1:J],1,K),spar_arr,N,J);
        Phi = sqrt(1/M)*randn(M,N);
        Y = Phi*S;
        
        tic
        S1 = SOMP(Y, Phi, K);
        toc
    
        timesaver1 = [timesaver1,toc];
        rmsesaver1 = [rmsesaver1,RMSE(S1,S)];

        tic
        S2 = MFOCUSS(Phi,Y,1e-10);
        toc

        timesaver2 = [timesaver2,toc];
        rmsesaver2 = [rmsesaver2,RMSE(S2,S)];
        
        tic
        S3 = spg_mmv(Phi,Y,0,options);
        toc

        timesaver3 = [timesaver3,toc];
        rmsesaver3 = [rmsesaver3,RMSE(S3,S)];

        tic
        S4 = sniht(Y,Phi,K);
        toc
    
        timesaver4 = [timesaver4,toc];
        rmsesaver4 = [rmsesaver4,RMSE(S4,S)];
        
        tic
        S5 = ADMML21_1(Y,Phi,1e-6,1e-5);
        toc

        timesaver5 = [timesaver5,toc];
        rmsesaver5 = [rmsesaver5,RMSE(S5,S)];
    
        tic
        S6 = MMV_ADMM_L20(Y,Phi,K,1);
        toc

        timesaver6 = [timesaver6,toc];
        rmsesaver6 = [rmsesaver6,RMSE(S6,S)];
    end

    Timesaver1 = [Timesaver1;timesaver1];
    Timesaver2 = [Timesaver2;timesaver2];
    Timesaver3 = [Timesaver3;timesaver3];
    Timesaver4 = [Timesaver4;timesaver4];
    Timesaver5 = [Timesaver5;timesaver5];
    Timesaver6 = [Timesaver6;timesaver6];

    RMSEsaver1 = [RMSEsaver1;rmsesaver1];
    RMSEsaver2 = [RMSEsaver2;rmsesaver2];
    RMSEsaver3 = [RMSEsaver3;rmsesaver3];
    RMSEsaver4 = [RMSEsaver4;rmsesaver4];
    RMSEsaver5 = [RMSEsaver5;rmsesaver5];
    RMSEsaver6 = [RMSEsaver6;rmsesaver6];
end

Avetime1 = mean(Timesaver1,2);
Avetime2 = mean(Timesaver2,2);
Avetime3 = mean(Timesaver3,2);
Avetime4 = mean(Timesaver4,2);
Avetime5 = mean(Timesaver5,2);
Avetime6 = mean(Timesaver6,2);

%plot
S = ['-ks';'-ko';'-kd';'-kv';'-k*';'-k+'];
figure;
ind = 1:20;
plot(ind,Avetime1,S(1,:),'color','b','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,Avetime2,S(2,:),'color','r','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,Avetime3,S(3,:),'color','g','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,Avetime4,S(4,:),'color','c','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,Avetime5,S(5,:),'color','m','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,Avetime6,S(6,:),'color','k','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold off;

xlim([1 20]);
y=cell(1,length(ind));
y(1:2:end)=num2cell(1:2:20);
set(gca,'xtick',ind,'xticklabel',y);
set(gca,'FontName','Times New Roman','Fontsize',11);
h=legend('SOMP','MFOCUSS','MMV-SPG','SNIHT','MMV-ADMM-L_{2,1}','MMV-ADMM-L_{2,0}','Location','north');
set(h,'FontSize',8);
xlabel('J');
ylabel('Time(s)');
%title('Time of different algorithms with different sensor number');

%successful percentage
epsilon = 1e-5;

SuccessPercentage1 = zeros(1,length(J_set));
SuccessPercentage2 = zeros(1,length(J_set));
SuccessPercentage3 = zeros(1,length(J_set));
SuccessPercentage4 = zeros(1,length(J_set));
SuccessPercentage5 = zeros(1,length(J_set));
SuccessPercentage6 = zeros(1,length(J_set));

for i = 1:length(J_set)
    flag1 = 0;
    flag2 = 0;
    flag3 = 0;
    flag4 = 0;
    flag5 = 0;
    flag6 = 0;
    for j = 1:CNT
        if RMSEsaver1(i,j)<epsilon
            flag1 = flag1 + 1;
        end
        if RMSEsaver2(i,j)<epsilon
            flag2 = flag2 + 1;
        end
        if RMSEsaver3(i,j)<epsilon
            flag3 = flag3 + 1;
        end
        if RMSEsaver4(i,j)<epsilon
            flag4 = flag4 + 1;
        end
        if RMSEsaver5(i,j)<epsilon
            flag5 = flag5 + 1;
        end
        if RMSEsaver6(i,j)<epsilon
            flag6 = flag6 + 1;
        end
    end

    SuccessPercentage1(i) = flag1/100;
    SuccessPercentage2(i) = flag2/100;
    SuccessPercentage3(i) = flag3/100;
    SuccessPercentage4(i) = flag4/100;
    SuccessPercentage5(i) = flag5/100;
    SuccessPercentage6(i) = flag6/100;
end

figure;
plot(ind,SuccessPercentage1,S(1,:),'color','b','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,SuccessPercentage2,S(2,:),'color','r','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,SuccessPercentage3,S(3,:),'color','g','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,SuccessPercentage4,S(4,:),'color','c','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,SuccessPercentage5,S(5,:),'color','m','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold on;
plot(ind,SuccessPercentage6,S(6,:),'color','k','Linewidth',1);
set(gca,'XTick',ind);
set(gca,'XTickLabel',J_set);
hold off;

xlim([1 20]);
y=cell(1,length(ind));
y(1:2:end)=num2cell(1:2:20);
set(gca,'xtick',ind,'xticklabel',y);
set(gca,'FontName','Times New Roman','Fontsize',11);
h=legend('SOMP','MFOCUSS','MMV-SPG','SNIHT','MMV-ADMM-L_{2,1}','MMV-ADMM-L_{2,0}','Location','southeast');
set(h,'FontSize',8);
xlabel('J');
ylabel('Percentage of successful recovery');
%title('Percentage of successful recovery of different algorithms with different sensor number');

save test5
load test5
