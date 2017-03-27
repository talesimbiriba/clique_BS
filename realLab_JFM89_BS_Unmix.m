%
%
%   Author: Tales Imbiriba   
%   Last Modified in November 2016


clear all;
close all; 


% loading data!
%load jfm89_oliv_enst_mix.mat
%load jfm89_oliv_mag_mix.mat
%load jfm89_oliv_arn_mix.mat
% load jfm89_arn_ens_mix.mat
load jfm89_oliv_arn_ens_mix.mat



[L,R] = size(M);

% FCLS
[~, a_fcls] = estimateLinearModel(y, M, true);


% GKKM
kbwkkm = 0.1006;
lambda = 2;
tic
% band selection
%[kkmBS] = kernelKMeansBandSelection(M, Nb_kkm, kbw);
[kkmBS] = kernelKMeansBandSelectionAIC(M,kbwkkm,lambda);
yr=y(kkmBS,:);
Mr=M(kkmBS,:);

% unmix!
%[a_kkm,~,r_kkm] = tskHype(yr, Mr,[],[],kbwkkm); 
[a_kkm,~,r_kkm] = tskHype_reducedData(yr, Mr,[],[],kbwkkm,[],M,kkmBS);

% time to BS + Unmix
kkmTime = toc;

%computing dictionary mu
Kg = computeKernelMatrix(Mr,kbwkkm);
mu_kkm = max(max(Kg-eye(size(Kg))));    
%[rmse_kkm, std_kkm] = RMSEAndSTDForMatrix(r_kkm,y);

%[rmse, stdd] = RMSEAndSTDForMatrix(a_kkm, a_mix);
[rmse, stdd] = RMSEAndSTDForMatrix(r_kkm, y);
fprintf('GKKM & %.4f $\\pm$ %.4f & %d & %d & %.4f\\\\ \n', rmse, stdd, kkmTime, length(kkmBS), mu_kkm);




% Fullband SK-Hype
kbw_skp = 0.00885444926741354;
tic
[a_skp,~,r_skp] = tskHype(y, M,[],[],kbw_skp); 
skpTime = toc;

disp('Results with Image reconstruction error')
fprintf('Strategy & RMSE $\\pm$ STD &  Time & $N_b$ & $\\mu$\\\\ \\hline\n');

%[rmse_skp, std_skp] = RMSEAndSTDForMatrix(a_skp,a_mix);
[rmse_skp, std_skp] = RMSEAndSTDForMatrix(r_skp,y);

fprintf('SK-Hype & %2.4f $\\pm$ %2.4f &  %2.4f & %d & -\\\\ \\hline\n', rmse_skp, std_skp, skpTime,L);




ms=[5 10 20 30];
%ms = [30];

% find Gaussian kernel bandwidth!
K_s1 = computeKernelMatrix(M,1);
c=0;
count =1;
for i=1:L-1
    for j=i+1:L
        c(count) = K_s1(i,j);
        count = count + 1;
    end
end


cliqueCBSTime = zeros(size(ms));
greedyCBSTime = zeros(size(ms));

mu_clique = zeros(size(ms));
mu_greedy = zeros(size(ms));

Nb_clique = zeros(size(ms));
Nb_greedy = zeros(size(ms));


i = 1;
%Opt_opt = optimset('Algorithm','interior-point','TolFun',1e-10);
Opt_opt = optimset('Algorithm','interior-point');

for m=ms    
    % m=10;   % number of desired bands
    mu_0 = 1/(m-1);
    [kbw,fval] = fmincon(@(kbw)(abs(mean(c.^(1/(kbw^2)))-mu_0)),1,[],[],[],[],1e-10,1e100,[],Opt_opt);

    KM = computeKernelMatrix(M,kbw);

    % clique (CCBS)
    tic
    % band selection
    [cliqueCBS] = clique_coherence_bandselection( KM, mu_0, [], 1 );
    yr=y(cliqueCBS,:);
    Mr=M(cliqueCBS,:);
    
    % unmix!
    %[a_clique,~,r_clique] = tskHype(yr, Mr,[],[],kbw); 
    [a_clique,~,r_clique] = tskHype_reducedData(yr, Mr,[],[],kbw,[],M,cliqueCBS); 
    
    % time to BS + Unmix
    cliqueCBSTime(i) = toc;
    
    Nb_clique = length(cliqueCBS);
    Nb_kkm = Nb_clique;
    
    %computing dictionary mu
    Kg = computeKernelMatrix(Mr,kbw);
    mu_clique(i) = max(max(Kg-eye(size(Kg))));
    %[rmse_clique, std_clique] = RMSEAndSTDForMatrix(r_clique,y);
    
    % greedy (GCBS)

    tic
    % band selection
    [greedyCBS] = buildDictionaryUsingCoherenceFactorKM(KM, mu_0);
    yr=y(greedyCBS,:);
    Mr=M(greedyCBS,:);
    
    % unmix!
    %[a_greedy,~,r_greedy]  = tskHype(yr, Mr,[],[],kbw); 
    [a_greedy,~,r_greedy]  = tskHype_reducedData(yr, Mr,[],[],kbw,[],M,greedyCBS); 
    
    % time to BS + Unmix
    greedyCBSTime(i) = toc;
    
    Nb_greedy = length(greedyCBS);
    
    %computing dictionary mu
    Kg = computeKernelMatrix(Mr,kbw);
    mu_greedy(i) = max(max(Kg-eye(size(Kg))));
    %[rmse_greedy, std_greedy] = RMSEAndSTDForMatrix(r_greedy,y);
    
    
    
    
    

    fprintf('$M = %d$, $\\mu_0 = %.4f$, $\\sigma = %.4f$\\\\ \\hline\n', m, mu_0, kbw);
    
%     printf('CCBS & %.4f $\\pm$ %.4f & %.2f & %d %.4f\n', rmse_clique, std_clique, Nb_clique, mu_clique);
%     printf('GCBS & %.4f $\\pm$ %.4f & %.2f & %d %.4f\n', rmse_greedy, std_greedy, Nb_greedy, mu_greedy);
%     printf('GKKM & %.4f $\\pm$ %.4f & %.2f & %d %.4f\n', rmse_kkm, std_kkm, Nb_kkm, mu_kkm);
    
    
    
    %[rmse, stdd] = RMSEAndSTDForMatrix(a_clique, a_mix);
    [rmse, stdd] = RMSEAndSTDForMatrix(r_clique, y);
    fprintf('CCBS & %.4f $\\pm$ %.4f & %.4f & %d & %.4f\\\\ \n', rmse, stdd, cliqueCBSTime(i), Nb_clique, mu_clique(i));
    %[rmse, stdd] = RMSEAndSTDForMatrix(a_greedy, a_mix);
    [rmse, stdd] = RMSEAndSTDForMatrix(r_greedy, y);
    fprintf('GCBS & %.4f $\\pm$ %.4f & %.4f & %d & %.4f\\\\ \n', rmse, stdd, greedyCBSTime(i),Nb_greedy, mu_greedy(i));
       
    
    i = i + 1;
    
    
    
%     figure
%     subplot(3,1,1)
%     plot(wlength, abs(mean(y'-r_skp')))
%     xlim([0.3 2.6])
%     ylim([0 0.03])
%     set(gca,'FontSize',14)
%     legend('SK-Hype','location', 'NorthEast')
%     
%     subplot(3,1,2)
%     plot(wlength, abs(mean(y'-r_kkm')))
%     xlim([0.3 2.6])
%     ylim([0 0.03])
%     set(gca,'FontSize',14)
%     legend('GKKM','location', 'NorthEast')
%     
%     subplot(3,1,3)
%     plot(wlength, abs(mean(y'-r_clique')))
%     xlim([0.3 2.6])
%     ylim([0 0.03])
%     set(gca,'FontSize',14)
%     legend('CCBS','location', 'NorthEast')
    
end


 
