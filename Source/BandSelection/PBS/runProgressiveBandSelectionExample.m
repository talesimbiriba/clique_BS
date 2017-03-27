

% clc;
clear all;
close all;


% verbosity: '0' print only the results '1' print all
verbosity = 0; 
R = 8;
% numOfPixels = 2000;
numOfPixels = 20;
%nruns = 100;
nruns = 1;
SNR = 21;


% M_database = 0;     %Pavia University data:
M_database = 3;    % Cedric/Jie Data (Cuprite Minerals)

%mixture model (1: linear, 2: bilinear, 3, pnmm)
model = 3;

% data generation 
[y,M,a,std_noise] = generate_image(M_database, R, model, numOfPixels, SNR);

L = size(M,1);


% print simulation config
disp(['Verbosity = ',num2str(verbosity)]);
if M_database==0
   disp('Using Pavia U. Endmembers') 
elseif M_database==3
    disp('Using Cuprite Endmembers') 
end
if model==1
    disp('Using the LMM model')
elseif model==2
    disp('Using the GBM model')
elseif model==3
    disp('Using the PNMM model')
end

disp(['R = ',num2str(R)]);
disp(['Runs = ',num2str(nruns)]);
disp(['NumberOfPixels = ',num2str(numOfPixels)]);
disp(['SNR = ',num2str(SNR)]);


%% Progressive BS

kbw_pbs = 0.1006 * 0.25;

metric = 0;     % variance
%metric = 1;     % entropy
%nb = 20;
nbs = [10 16 21 42];



for nb = nbs,
    pbsNb = 0;
    mu_pbs = 0;
    for i=1:nruns,
        tic
        [ pbsBS{i} ] = bandPrioritizationBS(y, metric, nb);
        pbsBSTime(i) = toc;

         yr=y(pbsBS{i},:);
         Mr=M(pbsBS{i},:);

        tic
        a_pbs = tskHype(yr, Mr,[],[],kbw_pbs); 
        skhype_pbs_time(i) = toc;


        [rmse_pbs(i),std_pbs(i)] = RMSEAndSTDForMatrix(a,a_pbs);

        Kg = computeKernelMatrix(Mr,kbw_pbs);
        mu_pbs = mu_pbs + max(max(Kg-eye(size(Kg))));
        pbsNb = pbsNb + length(pbsBS{i});

    end

    mu_pbs = mu_pbs/nruns;
    pbsNb=pbsNb/nruns;


    rmses = rmse_pbs; 
    stds = std_pbs;
    times = pbsBSTime + skhype_pbs_time;
    fprintf('PBS & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times), pbsNb, mu_pbs);                           
end
%         %% KKM - Kernel K-Means 
% 
% 
% if verbosity==1,
%     fprintf('Test5: Kernel K-Means \n')
% end
% 
% kbws2 = [kbw/2, kbw, 2*kbw, 10*kbw, 20*kbw];
% lambdas = [2 4 6];
% [kbw_kkm, lambda ] = searchForBestParametersKKM(y,M,a, lambdas, kbws2)
% 
% 
% mu_kkm=0;
% kkmNb=0;
% for i=1:nruns,
%     
%     tic            
%     [kkmBS{i}] = kernelKMeansBandSelectionAIC(M, kbw_kkm, lambda);
%     kkmBSTime(i) = toc;
%     
%     yr=y(kkmBS{i},:);
%     Mr=M(kkmBS{i},:);
% 
%     tic
%     a_kkm = tskHype(yr, Mr,[],[],kbw_kkm); 
%     skhype_kkm_time(i) = toc;
%     
%     
%     [rmse_kkm(i),std_kkm(i)] = RMSEAndSTDForMatrix(a,a_kkm);
% 
%     Kg = computeKernelMatrix(Mr,kbw_kkm);
%     mu_kkm = mu_kkm + max(max(Kg-eye(size(Kg))));
%     kkmNb = kkmNb + length(kkmBS{i});
% end
% 
% 
% 
% mu_kkm = mu_kkm/nruns;
% kkmNb=kkmNb/nruns;
% 
% 
% rmses = rmse_kkm; 
% stds = std_kkm;
% times = kkmBSTime + skhype_kkm_time;
% fprintf('GKKM & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times), kkmNb, mu_kkm);                           
% 
%     
% %% SK-Hype Simulation
%     
% %kbws2 = [kbw/2, kbw, 2*kbw, 10*kbw, 20*kbw];
% % the best results were obtained with:
% kbws2 = [kbw/2];
% 
% 
% for kbw_skp=kbws2    
%     for i=1:nruns,    
%         tic 
%     	a_skp= tskHype(y, M,[],[],kbw_skp);         
%         skhype_time(i) = toc;
%         
%         [rmse_skp(i),std_skp(i)] = RMSEAndSTDForMatrix(a,a_skp);
%     end
% 
%     rmses = rmse_skp;
%     stds = std_skp;
%     times = skhype_time;
%     fprintf('SK-Hype & %2.4f $\\pm$ %2.4f &  %2.4f $\\pm$ %2.4f & %d\n', mean(rmses), mean(stds), mean(times), std(times),L);
% 
% 
% end  % end for over kbws

% if model==3
%    save ~/compMethods_simCupPNMM.mat
% else
%     save ~/compMethods_simCupGBM.mat
% end

