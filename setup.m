% startup script 
%
% Copyright (c) by Tales Imbiriba 2015.

disp(['executing startup script...']);

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located

addpath(mydir(1:end-1))
addpath(genpath([mydir,'Source']))

stringCommand = ['cd(''',mydir,'Source/MaxCliqueExecutables/'')'];
fid=fopen('goToCliqueDir.m','wt');
fprintf(fid, '%s',stringCommand);
fclose(fid);
