function emebsddi_wrapper_fun( data )
%EMEBSDDI_WRAPPER_FUN
% Run EMsoft EMEBSDDI program to perform dictionary indexing
% Unlike EMsoft itself, this code can handle hexagonal grid data
% Original: 2/22/20 (Edward Pang, MIT)
% Change log:
% -4/20/21 ELP: fix .nml file so it works for EMsoft 4.0
%
%%% Input: 'data', a struct containing the following fields:
% %%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data.path = 'testdata2/';  % path within EMdatapathname to find .data/.ang files and export final data
% data.angfile_in = 'Scan1.ang';  % name of Hough .ang file (located in path)
% data.phaseid = 1;   % which phase to consider (as numbered in .ang file header)
%
% data.img_exp = 'Scan1.data';     % name of .data file containing all patterns (located in path)
% data.binning = 8;        % binning mode (1, 2, 4, 8)
% 
% data.outputname = 'test_RunEMEBSDDI_refine_Scan1';     % name of outputted files will be [outputname].h5/.ang/.nml
% 
% % Pattern center
% data.L = 2510.68;  % in microns
% data.xpc = 9.7354;    % in px
% data.ypc = 43.0779;   % in px
% 
% % Define orientation resolution (cubochoric N parameter)
% %   Larger value gives better orientation resolution but takes longer time to calculate.
% %   N=90: 1.5deg resolution, N=109: 1.25deg, N=137: 1deg. N>137 takes a long time (days) on most computers.
% data.N = 90;
% 
% % Other EMsoft parameters
% data.energymin = 15.0;   % energy range in the intensity summation (keV), must match master pattern
% data.energymax = 30.0;
% data.masterfile = '200106_Cu_fcc_30kV_60deg_master/Cu_fcc-master-30kV.h5';   % master pattern input file, path relative to EMdatapathname
%
% 
%%% Parameters you don't need to change often %%%
% % Paths for this computer
% data.homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
%
% % Detector parameters
% data.thetac = 8;   % tilt angle of the camera (positive below horizontal, degrees)
% data.delta = 7.4*480/416;   % CCD pixel size on the scintillator surface (microns)
% data.numsx = 416;    % number of CCD pixels along x and y
% data.numsy = 416;
% data.omega = 0;      % angle between normal of sample and detector
% data.maskpattern = 'y';  % apply circular mask to data? (y/n)
% data.maskradius = 104;   % RADIUS in pixels, AFTER binning operation
% 
% % Parameters for EMsoft EMEBSDDI
% data.scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
% data.gammavalue = 0.33;  % gamma correction factor
% data.nthreads = 8;
% data.hipassw = 0.05;        % hi pass filter w param; 0.05 is reasonable
% data.nregions = 10;      % # of regions for adaptive histogram equalization
% data.numdictsingle = 1024;   % number of dictionary files arranged in column for dp on GPU (multiples of 16 perform better)
% data.numexptsingle = 1024;   % number of experiment files "
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Extract input parameters from struct 'data'
path = data.path;
angfile_in = data.angfile_in;
img_exp = data.img_exp;
phaseid = data.phaseid;
outputname = data.outputname;
L = data.L;
xpc = data.xpc;
ypc = data.ypc;
N = data.N;
energymin = data.energymin;
energymax = data.energymax;
masterfile = data.masterfile;
binning = data.binning;
maskpattern = data.maskpattern;
maskradius = data.maskradius;
thetac = data.thetac;
delta = data.delta;
numsx = data.numsx;
numsy = data.numsy;
omega = data.omega;
homepath = data.homepath;
scalingmode = data.scalingmode;
gammavalue = data.gammavalue;
nthreads = data.nthreads;
hipassw = data.hipassw;
nregions = data.nregions;
numdictsingle = data.numdictsingle;
numexptsingle = data.numexptsingle;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


totaltime = tic;


% Name paths
angpath_in = fullfile(homepath,path,angfile_in);
angpath_out_short = fullfile(path,[outputname '.ang']);
angpath_out = fullfile(homepath,path,[outputname '.ang']);
nmlpath = fullfile(homepath,path,[outputname '.nml']);
h5path = fullfile(homepath,path,[outputname '.h5']);
h5pathshort = fullfile(path,[outputname '.h5']);
imgpathshort = fullfile(path,img_exp);


% Error checking
if exist(h5path,'file')==2
    error('File already exists! Rename .h5 output file.');
end
if exist(angpath_out,'file')==2
    error('File already exists! Rename .ang output file.');
end
if exist(nmlpath,'file')==2
    error('File already exists! Rename .nml output file.');
end


% read in .ang file
[~, x, y, IQ, ~, ~, ~, phaseinfo, grid, PChough] = loadang(angpath_in, 0);

% get phase info for phase of interest
for ii=1:size(phaseinfo,1)
    if phaseinfo{ii,1}==phaseid
        % store phase info
        phaseinfo_temp{1,1} = phaseinfo{ii,1};
        phaseinfo_temp{1,2} = phaseinfo{ii,2};
        phaseinfo_temp{1,3} = phaseinfo{ii,3};
        phaseinfo_temp{1,4} = phaseinfo{ii,4};
        break   % found it, don't need to keep checking
    end
end
if ~exist('phaseinfo_temp','var')
    error('Cannot find phaseid=%g.',phaseid);
end

% grid info
nrows = grid{6};
ncols_odd = grid{4};


% Create EMEBSDDI.nml file
fid = fopen(nmlpath,'w');
fprintf(fid,' &EBSDIndexingdata\n');
fprintf(fid,' indexingmode = ''dynamic'',\n');
fprintf(fid,' Notify = ''off'',\n');
fprintf(fid,' ipf_ht = %g,\n',nrows);       % height of data set in pattern input file
fprintf(fid,' ipf_wd = %g,\n',ncols_odd);   % width "
fprintf(fid,' ROI = 0 0 0 0,\n');   % leave all at 0 for full field of view
fprintf(fid,' stepX = 1.0,\n');     % X and Y sampling step sizes
fprintf(fid,' stepY = 1.0,\n');
fprintf(fid,' nnk = 50,\n');        % # of top dot products to save
fprintf(fid,' nnav = 20,\n');       % # of top matches to use for orientation averaging (<nnk)
fprintf(fid,' nosm = 20,\n');       % # of top matches to use for Orientation Similarity Map computation (<nnk)
% fprintf(fid,' nism = 5,\n');        % # of top matches to use for Indexing Success Map computation (<nnk)   %%% Not in EMsoft 4.0
% fprintf(fid,' isangle = 1.5,\n');   % indexing success threshold angle (deg)   %%% Not in EMsoft 4.0
fprintf(fid,' maskfile = ''undefined'',\n');
fprintf(fid,' maskpattern = ''%s'',\n',maskpattern);
fprintf(fid,' maskradius = %.0f,\n',maskradius);
fprintf(fid,' hipassw = %g,\n',hipassw);
fprintf(fid,' nregions = %g,\n',nregions);

fprintf(fid,' ncubochoric = %g,\n',N);
fprintf(fid,' L = %g,\n',L);
fprintf(fid,' thetac = %g,\n',thetac);
fprintf(fid,' delta = %g,\n',delta);
fprintf(fid,' numsx = %g,\n',numsx);
fprintf(fid,' numsy = %g,\n',numsy);
fprintf(fid,' xpc = %g,\n',xpc);
fprintf(fid,' ypc = %g,\n',ypc);
fprintf(fid,' omega = %g,\n',omega);
fprintf(fid,' energymin = %g,\n',energymin);
fprintf(fid,' energymax = %g,\n',energymax);
fprintf(fid,' spatialaverage = ''n'',\n');
fprintf(fid,' beamcurrent = 150.0,\n');
fprintf(fid,' dwelltime = 100.0,\n');
fprintf(fid,' binning = %g,\n',binning);
fprintf(fid,' scalingmode = ''%s'',\n',scalingmode);
fprintf(fid,' gammavalue = %g,\n',gammavalue);

fprintf(fid,' exptfile = ''%s'',\n',imgpathshort);
fprintf(fid,' inputtype = ''Binary'',\n');
fprintf(fid,' HDFstrings = '''' '''' '''' '''' '''' '''' '''' '''' '''' '''',\n');

fprintf(fid,' tmpfile = ''EMEBSDDict_tmp.data'',\n');
fprintf(fid,' keeptmpfile = ''n'',\n');
fprintf(fid,' datafile = ''%s'',\n',h5pathshort);
fprintf(fid,' ctffile = ''undefined'',\n');
fprintf(fid,' avctffile = ''undefined'',\n');
fprintf(fid,' eulerfile = ''undefined'',\n');

fprintf(fid,' dictfile = ''undefined'',\n');

fprintf(fid,' masterfile = ''%s'',\n',masterfile);

fprintf(fid,' numdictsingle = %g,\n',numdictsingle);
fprintf(fid,' numexptsingle = %g,\n',numexptsingle);
fprintf(fid,' nthreads = %g,\n',nthreads);
fprintf(fid,' platid = 1,\n');
fprintf(fid,' devid = 1,\n');
fprintf(fid,' multidevid = 1 0 0 0 0 0 0 0,\n');
fprintf(fid,' usenumd = 1,\n');

fprintf(fid,' /\n');
fclose(fid);


% Run EMsoft EMEBSDDI
system(sprintf('EMEBSDDI %s',nmlpath));     % print output to screen live



% Read in data
% sometimes can't find h5path immediately after running EMEBSDDI, these lines of code try to remedy that
t = tic;    
while exist(h5path,'file')==0 && toc(t)<60
    pause(0.1);     % pause 0.1 seconds and try again
end
if exist(h5path,'file')==2
    TopDotProductList = h5read(h5path,'/Scan 1/EBSD/Data/TopDotProductList');
    NumExptPatterns = h5read(h5path,'/Scan 1/EBSD/Data/NumExptPatterns');
    phi1 = h5read(h5path,'/Scan 1/EBSD/Data/Phi1');     % in rad
    PHI = h5read(h5path,'/Scan 1/EBSD/Data/Phi');     % in rad
    phi2 = h5read(h5path,'/Scan 1/EBSD/Data/Phi2');     % in rad
else
    error('h5path not found!');
end


% extract indexed solutions
dpbest = TopDotProductList(1,1:NumExptPatterns);
eulerbest = [phi1'; PHI'; phi2'];

    
% if hex grid, remove dummy data
if strcmp(grid{1},'HexGrid')
    % calculate indices of dummy data
    ix = 2*ncols_odd;
    dummy = ix:ix:NumExptPatterns;  % indices of dummy data

    % remove dummy data
    dpbest(dummy) = [];
    eulerbest(:,dummy) = [];
end



% calculate some things
fit = 1-dpbest';     % measure of OIM's fit for DI
CI = ones(NumExptPatterns,1);
phase = zeros(NumExptPatterns,1);   % OIM labels phase=0 in single phase .ang files
xstar = (xpc/numsx) + 0.5;
ystar = (ypc/numsy) + 0.5;
zstar = L/(numsx*delta);


%%% save data to struct for passing to makeang.m
data.eulerbest = eulerbest;
data.CI = CI;
data.fit = fit;
data.info = 'Indexed by EMsoft EMEBSDDI';
data.PC = [xstar ystar zstar PChough(4)];   % DI pattern center and WD (read from microscope ang)

% copy some things from Hough
data.IQ = IQ;
data.phase = phase;     % should be all zeros
data.x = x;     
data.y = y;
data.phaseinfo = phaseinfo_temp;
data.grid = grid;



%%% save .ang file
makeang(angpath_out, data, 1)
fprintf('Data stored in ang file: %s\n\n',angpath_out_short);



toc(totaltime)


