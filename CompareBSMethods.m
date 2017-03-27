

% clc;
clear all;
close all;


% verbosity: '0' print only the results '1' print all
verbosity = 0; 
R = 8;
% numOfPixels = 2000;
numOfPixels = 1000;
% nruns = 100;
nruns = 10;
SNR = 21;


% M_database = 0;     %Pavia University data:
M_database = 3;    % Cedric/Jie Data (Cuprite Minerals)

%mixture model (1: linear, 2: bilinear, 3, pnmm)
model = 3;

% data generation 
[y,M,a,std_noise]=generate_image(M_database, R, model, numOfPixels, SNR);

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



% c vector used to find kernel bandwidth!
K_s1 = computeKernelMatrix(M,1);
c=0;
count =1;
for i=1:L-1
    for j=i+1:L
        c(count) = K_s1(i,j);
        count = count + 1;
    end
end


% nBands = [5 10 20 30];
 nBands = [20];

for m=nBands,
% for m=20,

    % find Gaussian kernel bandwidth!
    %Opt_opt = optimoptions('fmincon','Algorithm','interior-point');
    Opt_opt = optimset('Algorithm','interior-point');
    % m=10;   % number of desired bands
    mu_0 = 1/(m-1);
    %[kbw,fval] = fmincon(@(kbw)(abs(mean(c.^(1/(kbw^2)))-mu_0)),0.4/m,[],[],[],[],1e-10,1e100,[],Opt_opt);
    [kbw,fval] = fmincon(@(kbw)(abs(mean(c.^(1/(kbw^2)))-mu_0)),1,[],[],[],[],1e-10,1e100,[],Opt_opt);


%     fprintf('kbw = %1.2f\n',kbw)
%     fprintf('m = %1.2f\n',m)
%     fprintf('mu_0 = %1.2f\n',mu_0)
%     fprintf('\n\nkbw = %1.2f\n',kbw);

    KM = computeKernelMatrix(M,kbw);
    MM = M*M';


    if verbosity==1,
        fprintf('Test1: several runs of the CCBS for a fixed M\n')
    end
    
    mu_clique = 0;
    for i=1:nruns,
        [cliqueCBS{i}, cliqueCBSTime(i)] = clique_coherence_bandselection( KM, mu_0, [], 1 );
       
        yr=y(cliqueCBS{i},:);
        Mr=M(cliqueCBS{i},:);

        % Unmix! 
        tic
        a_ccbs_f = tskHype(yr, Mr,[],[],kbw); 
        skhype_ccbs_f_time(i) = toc;

        % compute RMSE
        [rmse_ccbs_f(i),std_ccbs_f(i)] = RMSEAndSTDForMatrix(a,a_ccbs_f);
        
        Kg = computeKernelMatrix(Mr,kbw);
        mu_clique = mu_clique + max(max(Kg-eye(size(Kg))));
    end
    
    mu_clique = mu_clique/nruns;
    
    if verbosity==1,
        if sum(sum(repmat(cliqueCBS{1},1,nruns) - cell2mat(cliqueCBS))) == 0
            fprintf('All cliques are identical with clique size = %d.\n',length(cliqueCBS{1}))
        end
    end

    cliqueCBS_Nb = length(cliqueCBS{1});

    
    if verbosity==1
        rmses = rmse_ccbs_f;
        stds = std_ccbs_f;
        times = cliqueCBSTime + skhype_ccbs_f_time;

        fprintf('BS Strategy & RMSE $\\pm$ STD &  Time (Average $\\pm$ STD) & $N_b$ & $\\mu$\n ');
        fprintf('CCBS & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times),cliqueCBS_Nb, mu_clique);
        fprintf('\n')
    end

    
    %% sim 2 randomizing the rows of M
    if verbosity==1,
        fprintf('Test2: several runs of the CCBS for randomly rearrenged M\n')
    end

    %     matlabpool('close')
    %     matlabpool('open', nCores)
    cliqueCBS_random_Nb_vec = zeros(nruns,1);
    mu_clique_r = 0;
    for i=1:nruns,
        rss(:,i) = randsample(L,L);

        [cliqueCBS_random{i}, cliqueCBS_randomTime(i)] = clique_coherence_bandselection( KM(rss(:,i),rss(:,i)),mu_0,[],1);

        selBands = sort(rss(cliqueCBS_random{i},i));

        yr=y(selBands,:);
        Mr=M(selBands,:);

        % Unmix data!
        tic
        a_ccbs_r= tskHype(yr, Mr,[],[],kbw); 
        
        skhype_ccbs_r_time(i) = toc;

        [rmse_ccbs_r(i),std_ccbs_r(i)] = RMSEAndSTDForMatrix(a,a_ccbs_r);

        cliqueCBS_random_Nb_vec(i) = length(selBands);        
        
        Kg = computeKernelMatrix(Mr,kbw);
        mu_clique_r = mu_clique_r + max(max(Kg-eye(size(Kg))));
    end
    
    mu_clique_r = mu_clique_r/nruns;
    
    %cliqueCBS_random_Nb = size(cell2mat(cliqueCBS_random),1);
    cliqueCBS_random_Nb = mean(cliqueCBS_random_Nb_vec);
    cliqueCBS_random_Nb_std = std(cliqueCBS_random_Nb_vec);
    
    
    if verbosity==1
        rmses = rmse_ccbs_r;
        stds = std_ccbs_r;
        times = cliqueCBS_randomTime + skhype_ccbs_r_time;
        fprintf('CCBS (rand) & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times),cliqueCBS_random_Nb,mu_clique_r);
        fprintf('\n')
    end



    %% Greed Coherence BS (GCBS)
    if verbosity==1,
        fprintf('Test3: several runs of the GCBS for randomly rearrenged M\n')
    end

    mu_greed = 0;
    for i=1:nruns,
        tic 
        greedCBS{i} = buildDictionaryUsingCoherenceFactorKM( KM, mu_0);
        
        greedCBSTime(i) = toc;

        yr=y(greedCBS{i},:);
        Mr=M(greedCBS{i},:);
        
        tic
        a_gcbs_f = tskHype(yr, Mr,[],[],kbw); 
        skhype_gcbs_f_time(i) = toc;
        
        [rmse_gcbs_f(i),std_gcbs_f(i)] = RMSEAndSTDForMatrix(a,a_gcbs_f);
        
        Kg = computeKernelMatrix(Mr,kbw);
        mu_greed = mu_greed + max(max(Kg-eye(size(Kg))));
    end
    
    mu_greed = mu_greed/nruns;
    
    if verbosity==1,
        if sum(sum(repmat(greedCBS{1},1,nruns) - cell2mat(greedCBS)))==0
            fprintf('All set of selected bands are identical with Nb = %d.\n',size(cell2mat(greedCBS),1))
        end
    end

    greedCBS_Nb = size(cell2mat(greedCBS),1);

    
    if verbosity==1
        rmses = rmse_gcbs_f;
        stds = std_gcbs_f;
        times = greedCBSTime + skhype_gcbs_f_time;

        fprintf('GCBS & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times), greedCBS_Nb, mu_greed);
        fprintf('\n')
    end

    %%
    if verbosity==1,
        fprintf('Test4: several runs of the GCBS for randomly rearrenged M\n')
    end
    %     matlabpool('close')
    %     matlabpool('open', nCores)

    
    mu_greed_r=0;
    for i=1:nruns,
        rss(:,i) = randsample(L,L)';
        tic 
        greedCBS_random{i} = buildDictionaryUsingCoherenceFactorKM(KM(rss(:,i),rss(:,i)), mu_0);
        
        greedCBSTime_random(i) = toc;
        greedCBSNBands(i) = length(greedCBS_random{i});

        selBands = sort(rss(greedCBS_random{i},i));

        yr=y(selBands,:);
        Mr=M(selBands,:);
        
        tic
        a_gcbs_r= tskHype(yr, Mr,[],[],kbw); 
        skhype_gcbs_r_time(i) = toc;
        
        [rmse_gcbs_r(i),std_gcbs_r(i)] = RMSEAndSTDForMatrix(a,a_gcbs_r);

        Kg = computeKernelMatrix(Mr,kbw);
        mu_greed_r = mu_greed_r + max(max(Kg-eye(size(Kg))));
    end

    mu_greed_r = mu_greed_r/nruns;
    
    
    greedCBS_random_Nb = mean(greedCBSNBands);
    greedCBS_random_Nb_std = std(greedCBSNBands);
    if verbosity==1
        rmses = rmse_gcbs_r;
        stds = std_gcbs_r;
        times = greedCBSTime_random + skhype_gcbs_r_time;

        fprintf('Selected different set of bands with different Nb. E{Nb} = %2.2f, and std{Nb} = %2.2f \n',mean(greedCBSNBands),std(greedCBSNBands))
        fprintf('GCBS (rand) & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %2.2f $\\pm$ %2.2f & 2.4f\n', mean(rmses), mean(stds), mean(times), std(times), greedCBS_random_Nb, greedCBS_random_Nb_std, mu_greed_r);
        fprintf('\n')
    end


    %% Print Results 





    fprintf('\n\n Results for m = %d, mu_0 = %.4f, sigma = %.4f : \n\n',m,mu_0,kbw)
    fprintf('BS Strategy & RMSE $\\pm$ STD &  Time (Average $\\pm$ STD) & $N_b$ & $\\mu$\n');
    rmses = rmse_ccbs_f;
    stds = std_ccbs_f;
    times = cliqueCBSTime + skhype_ccbs_f_time;
    fprintf('CCBS & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times),cliqueCBS_Nb, mu_clique);

    rmses = rmse_ccbs_r;
    stds = std_ccbs_r;
    times = cliqueCBS_randomTime + skhype_ccbs_r_time;
    fprintf('CCBS (rand) & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f &  %2.2f $\\pm$ %2.2f & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times),cliqueCBS_random_Nb,cliqueCBS_random_Nb_std, mu_clique_r);

    rmses = rmse_gcbs_f;
    stds = std_gcbs_f;
    times = greedCBSTime + skhype_gcbs_f_time;
    fprintf('GCBS & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times), greedCBS_Nb, mu_greed );

    rmses = rmse_gcbs_r;
    stds = std_gcbs_r;
    times = greedCBSTime_random + skhype_gcbs_r_time;
    fprintf('GCBS (rand) & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %2.2f $\\pm$ %2.2f & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times), greedCBS_random_Nb, greedCBS_random_Nb_std, mu_greed_r);    
    


    fprintf('\n')



end
    

        %% KKM - Kernel K-Means 


if verbosity==1,
    fprintf('Test5: Kernel K-Means \n')
end

kbws2 = [kbw/2, kbw, 2*kbw, 10*kbw, 20*kbw];
lambdas = [2 4 6];
[kbw_kkm, lambda ] = searchForBestParametersKKM(y,M,a, lambdas, kbws2)


mu_kkm=0;
kkmNb=0;
for i=1:nruns,
    
    tic            
    [kkmBS{i}] = kernelKMeansBandSelectionAIC(M, kbw_kkm, lambda);
    kkmBSTime(i) = toc;
    
    yr=y(kkmBS{i},:);
    Mr=M(kkmBS{i},:);

    tic
    a_kkm = tskHype(yr, Mr,[],[],kbw_kkm); 
    skhype_kkm_time(i) = toc;
    
    
    [rmse_kkm(i),std_kkm(i)] = RMSEAndSTDForMatrix(a,a_kkm);

    Kg = computeKernelMatrix(Mr,kbw_kkm);
    mu_kkm = mu_kkm + max(max(Kg-eye(size(Kg))));
    kkmNb = kkmNb + length(kkmBS{i});
end



mu_kkm = mu_kkm/nruns;
kkmNb=kkmNb/nruns;


rmses = rmse_kkm; 
stds = std_kkm;
times = kkmBSTime + skhype_kkm_time;
fprintf('GKKM & %2.4f $\\pm$ %2.4f &  %2.2f $\\pm$ %2.2f & %d & %2.4f\n', mean(rmses), mean(stds), mean(times), std(times), kkmNb, mu_kkm);                           

    
%% SK-Hype Simulation
    
%kbws2 = [kbw/2, kbw, 2*kbw, 10*kbw, 20*kbw];
% the best results were obtained with:
kbws2 = [kbw/2];


for kbw_skp=kbws2    
    for i=1:nruns,    
        tic 
    	a_skp= tskHype(y, M,[],[],kbw_skp);         
        skhype_time(i) = toc;
        
        [rmse_skp(i),std_skp(i)] = RMSEAndSTDForMatrix(a,a_skp);
    end

    rmses = rmse_skp;
    stds = std_skp;
    times = skhype_time;
    fprintf('SK-Hype & %2.4f $\\pm$ %2.4f &  %2.4f $\\pm$ %2.4f & %d\n', mean(rmses), mean(stds), mean(times), std(times),L);


end  % end for over kbws

% if model==3
%    save ~/compMethods_simCupPNMM.mat
% else
%     save ~/compMethods_simCupGBM.mat
% end

