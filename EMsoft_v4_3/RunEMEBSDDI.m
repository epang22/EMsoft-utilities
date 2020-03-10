%RUNEMEBSDDI
% Perform dictionary indexing by calling EMsoft EMEBSDDI program
% Unlike EMsoft itself, this code can handle hexagonal grid data
% Fill in INPUT PARAMETERS section with desired parameters
% 2/22/20 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.path = 'testdata2/';  % path within EMdatapathname to find .data/.ang files and export final data
data.angfile_in = 'Scan1.ang';  % name of Hough .ang file (located in path)
data.phaseid = 1;   % which phase to consider (as numbered in .ang file header)

data.img_exp = 'Scan1.data';     % name of .data file containing all patterns (located in path)
data.binning = 8;        % binning mode (1, 2, 4, 8)

data.outputname = 'test_RunEMEBSDDI_Scan1';     % name of outputted files will be [outputname].h5/.ang/.nml

% Pattern center
data.L = 2510.68;  % in microns
data.xpc = 9.7354;    % in px
data.ypc = 43.0779;   % in px

% Define orientation resolution (cubochoric N parameter)
%   Larger value gives better orientation resolution but takes longer time to calculate.
%   N=90: 1.5deg resolution, N=109: 1.25deg, N=137: 1deg. N>137 takes a long time (days) on most computers.
data.N = 90;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %
        


%%% Run EMEBSDDI
emebsddi_wrapper_fun(data)


