
clear all
close all
clc;



% cuprite img used in detection paper Imbiriba et al. "Nonparametric
% Detection of Nonlinearly Mixed Pixels and Endmember Estimation in Hyperspectral Images"
imgNameAndPath =  'smallCupriteIMG.mat';

%load HyperCube smallIMG, 
% pixel Matrix Y
load(imgNameAndPath)

figure;
plotBands = [20 100 150];
imagesc(smallIMG(:,:,plotBands)./(max(max(max(smallIMG(:,:,plotBands ))))))

[L,N] = size(Y);

R = 5;
M = hyperVca(Y,R);


% y=Y(:,1:200);
y=Y;
%% Unmix! 

kbw_skp = 0.00885444926741354;
tic
[a_skp,~,r_skp] = tskHype(y, M,[],[],kbw_skp); 
skpTime = toc;

disp('Results with Image reconstruction error')
fprintf('Strategy & RMSE $\\pm$ STD &  Time & $N_b$ & $\\mu$\n');

[rmse_skp, std_skp] = RMSEAndSTDForMatrix(r_skp,y);

fprintf('SK-Hype & %2.4f $\\pm$ %2.4f &  %2.4f & %d & -\n', rmse_skp, std_skp, skpTime,L);




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
[a_kkm,~,r_kkm] = tskHype(yr, Mr,[],[],kbwkkm); 

% time to BS + Unmix
kkmTime = toc;

%computing dictionary mu
Kg = computeKernelMatrix(Mr,kbwkkm);
mu_kkm = max(max(Kg-eye(size(Kg))));    
%[rmse_kkm, std_kkm] = RMSEAndSTDForMatrix(r_kkm,y);

[rmse, stdd] = RMSEAndSTDForMatrix(a_kkm, a_skp);
fprintf('GKKM & %.4f $\\pm$ %.4f & %d & %d & %.4f\\\\ \n', rmse, stdd, kkmTime, length(kkmBS), mu_kkm);




% BS 
ms=[5 10 20 30];

% ms=[30];

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
for m=ms    
    Opt_opt = optimset('Algorithm','interior-point');
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
    [a_clique,~,r_clique] = tskHype(yr, Mr,[],[],kbw); 
    
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
    [a_greedy,~,r_greedy]  = tskHype(yr, Mr,[],[],kbw); 
    
    % time to BS + Unmix
    greedyCBSTime(i) = toc;
    
    Nb_greedy = length(greedyCBS);
    
    %computing dictionary mu
    Kg = computeKernelMatrix(Mr,kbw);
    mu_greedy(i) = max(max(Kg-eye(size(Kg))));
    %[rmse_greedy, std_greedy] = RMSEAndSTDForMatrix(r_greedy,y);
    
    
    
    
    

    fprintf('$m = %d$, $mu_0 = %.4f$, $\\sigma = %.4f$\n', m, mu_0, kbw);
    
%     printf('CCBS & %.4f $\\pm$ %.4f & %.2f & %d %.4f\n', rmse_clique, std_clique, Nb_clique, mu_clique);
%     printf('GCBS & %.4f $\\pm$ %.4f & %.2f & %d %.4f\n', rmse_greedy, std_greedy, Nb_greedy, mu_greedy);
%     printf('GKKM & %.4f $\\pm$ %.4f & %.2f & %d %.4f\n', rmse_kkm, std_kkm, Nb_kkm, mu_kkm);
    
    
    fprintf('\nCOMPARING ABUNDANCES BS AND SKHYPE\n\n')
    
    [rmse, stdd] = RMSEAndSTDForMatrix(a_clique, a_skp);
    fprintf('CCBS & %.4f $\\pm$ %.4f & %.4f & %d & %.4f\\\\ \n', rmse, stdd, cliqueCBSTime(i), Nb_clique, mu_clique(i));
    [rmse, stdd] = RMSEAndSTDForMatrix(a_greedy, a_skp);
    fprintf('GCBS & %.4f $\\pm$ %.4f & %.4f & %d & %.4f\\\\ \n', rmse, stdd, greedyCBSTime(i),Nb_greedy, mu_greedy(i));
       
    
    i = i + 1;
    
end












