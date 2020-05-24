%RUNEMEBSDDI_REFINE
% Perform dictionary indexing by calling EMsoft EMEBSDDI 
% program and then automatically starts a refinement run using
% EMFitOrientation
% Unlike EMsoft itself, this code can handle hexagonal grid data
% Can supply optional pseudosymmetry variants to check
% Fill in INPUT PARAMETERS section with desired parameters
% 2/23/20 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.path = 'testdata2/';  % path within EMdatapathname to find .data/.ang files and export final data
data.angfile_in = 'Scan1.ang';  % name of Hough .ang file (located in path)
data.phaseid = 1;   % which phase to consider (as numbered in .ang file header)

data.img_exp = 'Scan1.data';     % name of .data file containing all patterns (located in path)
data.binning = 8;        % binning mode (1, 2, 4, 8)

data.outputname = 'test_RunEMEBSDDI_refine_Scan1';
    % name of outputted intermediate files will be [outputname].h5/.ang/.ctf/.nml (located in path)
    % final output file will be [outputname]_EMFitOrientation.ang (located in path)

% Pattern center
data.L = 2510.68;  % in microns
data.xpc = 9.7354;    % in px
data.ypc = 43.0779;   % in px

% Define orientation resolution (cubochoric N parameter)
%   Larger value gives better orientation resolution but takes longer time to calculate.
%   N=90: 1.5deg resolution, N=109: 1.25deg, N=137: 1deg. N>137 takes a long time (days) on most computers.
data.N = 90;

% % Define pseudosymmetry variants to check (comment out if you don't want to check pseudosym)
% pseudosym = [
%     1 1 0 90;
%     1 -1 0 90
%     ];     % [axis_x axis_y axis_z angle(deg)]

% Other EMsoft parameters
data.energymin = 15.0;   % energy range in the intensity summation (keV), must match master pattern
data.energymax = 30.0;
data.masterfile = '200106_Cu_fcc_30kV_60deg_master/Cu_fcc-master-30kV.h5';   % master pattern input file, path relative to EMdatapathname



%%% Parameters you don't need to change often %%%
% Paths for this computer
data.homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)

% Detector parameters
data.thetac = 8;   % tilt angle of the camera (positive below horizontal, degrees)
data.delta = 7.4*480/416;   % CCD pixel size on the scintillator surface (microns)
data.numsx = 416;    % number of CCD pixels along x and y
data.numsy = 416;
data.omega = 0;      % angle between normal of sample and detector
data.maskpattern = 'y';  % apply circular mask to data? (y/n)
data.maskradius = min(data.numsx,data.numsy)/(2*data.binning);   % RADIUS in pixels, AFTER binning operation

% Parameters for EMsoft EMEBSDDI
data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
data.gammavalue = 0.33;  % gamma correction factor
data.nthreads = 8;
data.hipassw = 0.05;        % hi pass filter w param; 0.05 is reasonable
data.nregions = 10;      % # of regions for adaptive histogram equalization
data.numdictsingle = 1024;   % number of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
data.numexptsingle = 1024;   % number of experiment files "

% Parameters for EMsoft EMFitOrientation
step = 0.03;    %  max step size to take in homochoric space during the refinement (BOBYQA)
    % step = (.75*(angle-sin(angle)))^(1/3)
    % default: step=0.03, which corresponds to 4.16deg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



%%% error checking for EMFitOrientation
% Name paths
anginputpath = fullfile(data.homepath,data.path,[data.outputname '.ang']);
h5inputpath_short = fullfile(data.path,[data.outputname '.h5']);

ctfpath_short = fullfile(data.path,[data.outputname '_EMFitOrientation.ctf']);
ctfpath = fullfile(data.homepath,data.path,[data.outputname '_EMFitOrientation.ctf']);
angoutputpath = fullfile(data.homepath,data.path,[data.outputname '_EMFitOrientation.ang']);
nmlpath = fullfile(data.homepath,data.path,[data.outputname '_EMFitOrientation.nml']);

if exist('pseudosym','var')
    pspath_short = fullfile(data.path,[data.outputname '_PSvar.txt']);
    pspath = fullfile(data.homepath,data.path,[data.outputname '_PSvar.txt']);
end


% Error checking
if exist(ctfpath,'file')==2
    error('File already exists! Rename output file.');
end
if exist(angoutputpath,'file')==2
    error('File already exists! Rename output file.');
end
if exist(nmlpath,'file')==2
    error('File already exists! Rename output file.');
end
if exist('pseudosym','var')
    if exist(pspath,'file')==2
        error('File already exists! Rename output file.');
    end
end


%%% Run EMEBSDDI
emebsddi_wrapper_fun(data)



%%% EMFitOrientation
fprintf('\nStarting refinement run...\n');
tic


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
fprintf(fid,' nthreads = %g,\n',data.nthreads);
fprintf(fid,' dotproductfile = ''%s'',\n',h5inputpath_short);
fprintf(fid,' ctffile = ''%s'',\n',ctfpath_short);
fprintf(fid,' modality = ''EBSD'',\n');
fprintf(fid,' inRAM = .FALSE.,\n');
fprintf(fid,' matchdepth = 1,\n');
fprintf(fid,' method = ''FIT'',\n');
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

