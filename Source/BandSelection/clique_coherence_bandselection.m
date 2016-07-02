function [ selectedBands, bs_time ] = clique_coherence_bandselection( KM, nu_max, pathToDir, algorithm )
% [ selectedBands, bs_time ] = clique_coherence_bandselection( KM, nu_max,
% pathToDir, algorithm )
%   
%   'KM' is the kernel matrix 
%   'nu_max' is the coherence parameter 
%   'pathToDir' is the string with the path to the directory where the
%   executables are. 
%   'algorithm' is '1' (default) for the maxCLQ and '2' for the Cliquer

if nargin < 4
    algorithm = 1;
end
    
if nargin <3
    pathToDir = [];
end



% if isempty(pathToDir)
%     if algorithm==1
%         pathToDir = '~/Dropbox/matlab/HyperspectralCode/MaxCliqueExecutables/';
%     else
%     	pathToDir = '~/Work/matlab/temp/cliquer-1.21/';
%     end
% end
% 
% 
% cd(pathToDir)

goToCliqueDir;


%system('rm k_graph.a k_bs k_bs2');
fid = fopen('k_graph.a','wt');

L = size(KM,1);
fprintf(fid, 'p edge %d %d\n', L, length(find(KM < nu_max)));
for i=1:L
    for j=i+1:L
        if (KM(i,j)<nu_max)
            fprintf(fid, 'e %d %d\n', i, j);
        end
    end
end

fclose(fid);

if algorithm==1
    tic 
    system('./MaxCLQictai10Linux k_graph.a > k_bs');
    bs_time = toc;
    system('tail -n3 k_bs | head -n1  > k_bs2');
else
    tic 
    system('./cl -q -q --maximal k_graph.a > k_bs');
    bs_time = toc;
    system('cat k_bs | awk ''{$1=$2=""; print $0}''>k_bs2');
end

selectedBands = load('k_bs2');
selectedBands  = sort(selectedBands(:));
system('rm k_graph.a k_bs k_bs2 resultsOfiMCQ+MSZ5+TwoIset17');
end

