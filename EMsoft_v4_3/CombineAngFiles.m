%COMBINEANGFILES
% Combine data from multiple .ang files (one phase each) outputted by
% EBSDrefine or EMsoft-utilities, selecting the phase with the highest NDP
% Fill in INPUT PARAMETERS section with desired parameters
% 2/15/20 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List paths (within homepath) of each .ang file to read in
inputfile = {
    'testdata2/tetragonal.ang';
    'testdata2/monoclinic.ang'
    };

% Path (within homepath) to output combined .ang file
outputfile = 'testdata2/tetragonal+monoclinic.ang';


%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



% Error checking
angpath = fullfile(homepath,outputfile);    % full path to ang output file
if exist(angpath,'file')==2
    error('File already exists! Rename .ang output file.');
end


% read in ang files
Nfiles = length(inputfile);  % number of files being combined
for ii=1:Nfiles
    inputpath = fullfile(homepath,inputfile{ii});
    [euler, x, y, IQ, CI, ~, fit, phaseinfo, grid, PC] = loadang(inputpath,0);
    euler_all(:,3*ii-2:3*ii) = euler;
    x_all(:,ii) = x;
    y_all(:,ii) = y;
    IQ_all(:,ii) = IQ;
    CI_all(:,ii) = CI;
    fit_all(:,ii) = fit;
    phaseinfo_all{ii,1} = ii;
    phaseinfo_all{ii,2} = phaseinfo{2};
    phaseinfo_all{ii,3} = phaseinfo{3};
    phaseinfo_all{ii,4} = phaseinfo{4};
    grid_all{ii} = grid;
    PC_all(ii,:) = PC;
end

% check that data match up
gridinfomatch = 1;  % counter, change to zero if grid info doesn't match, to not repeat warnings
PCinfomatch = 1;  % counter, change to zero if PC info doesn't match, to not repeat warnings
for ii=2:Nfiles
    % check x,y positions
    if sum(x_all(:,1)~=x_all(:,ii)) > 0
        error('X values do not match.');
    end
    if sum(y_all(:,1)~=y_all(:,ii)) > 0
        error('Y values do not match.');
    end
    
    % check grid info
    if gridinfomatch==1     % check if up until now, grid info is same. if you know grid info is different, don't need to keep checking
        if ~strcmp(grid_all{1}{1},grid_all{ii}{1}) || grid_all{1}{2}~=grid_all{ii}{2} ||...
                grid_all{1}{3}~=grid_all{ii}{3} || grid_all{1}{4}~=grid_all{ii}{4} ||...
                grid_all{1}{5}~=grid_all{ii}{5} || grid_all{1}{6}~=grid_all{ii}{6}
            gridinfomatch = 0;  % change counter, now code knows grid info is different
            warning('Grids do not match. Using grid info from file 1: %s',inputfile{1});
        end
    end
    
    % check PC info
    if PCinfomatch==1     % check if up until now, PC info is same. if you know PC info is different, don't need to keep checking
        if sum(abs(PC_all(ii,:)-PC_all(1,:)))>0
            PCinfomatch = 0;  % change counter, now code knows PC info is different
            warning('PC does not match. Using PC info from file 1: %s',inputfile{1});
        end
    end
end


% Initialize vectors to store best solution and associated parameters
N = length(x);  % number of map points
CIbest = zeros(N,1);   
IQbest = zeros(N,1);
eulerbest = zeros(N,3);
fitbest = zeros(N,1);
phasebest = zeros(N,1);


% Determine best phase (lowest fit=highest dp) for each map pixel
for ii=1:N
    fit_temp = fit_all(ii,:);
    [m,index] = min(fit_temp);
    fitbest(ii) = m;
    phasebest(ii) = index;
    IQbest(ii) = IQ_all(ii,index);
    eulerbest(ii,:) = euler_all(ii,3*index-2:3*index);
    
    % compute CI
    fit_temp2 = fit_temp; fit_temp(index)=[];   % fits excluding best
    fitdiff = max(fit_temp2) - m;   % dp diff btw best and second best phase
    CIbest(ii) = min(CI_all(ii,index),fitdiff);     % CI is min btw single phase CI and fitdiff
end


% rearrange phaseinfo in order in which they appear
order = unique(phasebest,'stable');     % order in which phases appear
for ii=1:Nfiles
    phaseinfo_all_sort{ii,1} = phaseinfo_all{order(ii),1};
    phaseinfo_all_sort{ii,2} = phaseinfo_all{order(ii),2};
    phaseinfo_all_sort{ii,3} = phaseinfo_all{order(ii),3};
    phaseinfo_all_sort{ii,4} = phaseinfo_all{order(ii),4};
end


% save data to struct to pass to makeang.m
data.eulerbest = eulerbest'*pi/180;
data.CI = CIbest';
data.fit = fitbest';
data.IQ = IQbest';
data.info = 'Indexed by EBSDrefine (employing EMsoft EMEBSDDI)';
data.PC = PC_all(1,:);
data.phase = phasebest';
data.x = x';
data.y = y';
data.phaseinfo = phaseinfo_all_sort;
data.grid = grid_all{1};


% make ang file
makeang(angpath, data, 0);
fprintf('Data stored in ang file: %s\n',angpath);


