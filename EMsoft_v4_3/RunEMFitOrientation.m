%RUNEMFITORIENTATION
% Refine dictionary indexing data from EMsoft EMEBSDDI using EMsoft EMFitOrientation
% Unlike EMsoft itself, this code can handle hexagonal grid data
% Can supply optional pseudosymmetry variants to check
% Fill in INPUT PARAMETERS section with desired parameters
% 1/28/20 (Edward Pang, MIT)

clear
tic

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'testdata2/';  % path within EMdatapathname to find .data/.ang files and export final data

h5input = 'test_RunEMEBSDDI_Scan1.h5';  % name of .h5 file outputted by RunEMEBSDDI (located in path)
anginput = 'test_RunEMEBSDDI_Scan1.ang';    % name of .ang file outputted by RunEMEBSDDI (located in path)
    
outputname = 'test_RunEMEBSDDI_Scan1_EMFitOrientation';     % name of outputted files will be [outputname].ctf/.ang/.nml

% % Define pseudosymmetry variants to check (comment out if you don't want to check pseudosym)
% pseudosym = [
%     1 1 0 90;
%     1 -1 0 90
%     ];     % [axis_x axis_y axis_z angle(deg)]



%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)

% Parameters for EMsoft EMFitOrientation
nthreads = 8;
step = 0.03;    %  max step size to take in homochoric space during the refinement (BOBYQA)
    % step = (.75*(angle-sin(angle)))^(1/3)
    % default: step=0.03, which corresponds to 4.16deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %
        


% Name paths
anginputpath = fullfile(homepath,path,anginput);
h5inputpath_short = fullfile(path,h5input);

ctfpath_short = fullfile(path,[outputname '.ctf']);
ctfpath = fullfile(homepath,path,[outputname '.ctf']);
angoutputpath = fullfile(homepath,path,[outputname '.ang']);
nmlpath = fullfile(homepath,path,[outputname '.nml']);

if exist('pseudosym','var')
    pspath_short = fullfile(path,[outputname '_PSvar.txt']);
    pspath = fullfile(homepath,path,[outputname '_PSvar.txt']);
end


% Error checking
if exist(ctfpath,'file')==2
    error('File already exists! Rename .h5 output file.');
end
if exist(angoutputpath,'file')==2
    error('File already exists! Rename .ang output file.');
end
if exist(nmlpath,'file')==2
    error('File already exists! Rename .nml output file.');
end
if exist('pseudosym','var')
    if exist(pspath,'file')==2
        error('File already exists! Rename .txt output file.');
    end
end


% create PSvar file if specified
if exist('pseudosym','var')
    Npseudo = size(pseudosym,1);    % number of pseudosym vars to check
    fid = fopen(pspath,'w');
    fprintf(fid,'ax\n');
    fprintf(fid,'%g\n',Npseudo);
    for ii=1:Npseudo
        axis = pseudosym(ii,1:3);
        axis_unit = axis/norm(axis);    % make unit vector
        fprintf(fid,'%12.9f%12.9f%12.9f%12.9f\n',axis_unit(1),axis_unit(2),axis_unit(3),pseudosym(ii,4));
    end
    fclose(fid);
end


% Create EMFitOrientation.nml file
fid = fopen(nmlpath,'w');
fprintf(fid,' &RefineOrientations\n');
fprintf(fid,' nthreads = %g,\n',nthreads);
fprintf(fid,' dotproductfile = ''%s'',\n',h5inputpath_short);
fprintf(fid,' ctffile = ''%s'',\n',ctfpath_short);
fprintf(fid,' tmpfile = ''EMEBSDDict_tmp.data'',\n');
fprintf(fid,' modality = ''EBSD'',\n');
fprintf(fid,' inRAM = .FALSE.,\n');
fprintf(fid,' matchdepth = 1,\n');
fprintf(fid,' method = ''FIT'',\n');
fprintf(fid,' niter = 1,\n');
fprintf(fid,' nmis = 1,\n');
fprintf(fid,' step = %g,\n',step);
if exist('pseudosym','var')
    fprintf(fid,' PSvariantfile = ''%s''\n',pspath_short);
else
    fprintf(fid,' PSvariantfile = ''undefined''\n');
end
fprintf(fid,' /\n');
fclose(fid);


% Run EMFitOrientation
system(sprintf('EMFitOrientation %s',nmlpath));     % print output to screen live


% read in ang file data
[~,x,y,IQ,~,phase,~,phaseinfo,grid,PChough] = loadang(anginputpath,1);
N = size(x,1);     % number of map points


% read in ctf file data
[eulerbest, ~, ~, dp] = loadctf(ctfpath,1);


% figure out some grid info
ipf_ht = grid{6};
ipf_wd = grid{4};
if strcmp(grid{1},'HexGrid')
    issquare = 0;
else
    issquare = 1;   % square grid
end


% if data hexagonal, remove dummy data from ctf data
if issquare==0
    % figure out dummy indices
    idummy = [];
    for ii=2:2:ipf_ht
        itemp2 = 2*ipf_wd + 2*ipf_wd*(ii/2-1);
        idummy = [idummy itemp2];
    end
    
    % remove data
    eulerbest(idummy,:) = [];
    dp(idummy) = [];
end



% calculate some things
fit = 1-dp;     % measure of OIM's fit for DI
eulerbest = eulerbest*pi/180;   % convert to rad for storing
CI = ones(N,1);


%%% save data to struct for passing to makeang.m
data.eulerbest = eulerbest';
data.CI = CI';
data.fit = fit';
data.info = 'Refined by EMFitOrientation';

% copy some things from input ang
if ~isempty(PChough)
    data.PC = PChough;
end
data.IQ = IQ';
data.phase = phase';    % should all be zeros
data.x = x';
data.y = y';
data.phaseinfo = phaseinfo;
if ~isempty(grid)
    data.grid = grid;
end


% make ang file
makeang(angoutputpath, data, 1);
fprintf('Data stored in ang file: %s\n',angoutputpath);




toc


