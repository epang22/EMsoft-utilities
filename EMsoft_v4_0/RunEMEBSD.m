%RUNEMEBSD
% MATLAB wrapper for EMsoft's EMEBSD program for pattern simulation
% Run EMsoft EMEBSD program for a number of orientations, extract images, 
% and save them to disk
% Fill in INPUT PARAMETERS section with desired parameters
% 1/8/20 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'testdata2/';  % path within EMdatapathname to export images (create path before running)
outputname = 'test_RunEMEBSD';     % name of outputted image files will be [outputname]_#.png (# will be 1,2,3,etc.)

% Pattern center
L = 15000;  % in microns
xpc = 0;    % in px
ypc = 80;   % in px

% Define euler angles (each row phi1, PHI, phi2 in degrees)
euler = [
	0 0 0
	96.4666   27.6223  102.2867
    ];

% binning mode (1, 2, 4, 8)
binning = 2;

% EMsoft master pattern parameters
energymin = 15.0;   % energy range in the intensity summation (keV), must match master pattern
energymax = 30.0;
masterfile = '200106_Cu_fcc_30kV_60deg_master/Cu_fcc-master-30kV.h5';   % master pattern input file, path relative to EMdatapathname
      


%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)

% Detector parameters
thetac = 10;   % tilt angle of the camera (positive below horizontal, degrees)
delta = 59.2;   % CCD pixel size on the scintillator surface (microns)
numsx = 480;    % number of CCD pixels along x and y
numsy = 480;
omega = 0;      % angle between normal of sample and detector

% Pattern simulation parameters
imagetype = '.png';     % type of image file to save as
eulerconvention = 'tsl';    % 'tsl' (EDAX) or 'hkl' (Oxford) Euler angle convention parameter
scalingmode = 'gam'; % intensity scaling mode: 'not'=no scaling, 'lin'=linear, 'gam'=gamma
gammavalue = 0.33;  % gamma correction factor
nthreads = 1;
poisson = 'n';      % include poisson noise? (y/n)
maskpattern = 'n';  % apply circular mask to data? (y/n)
maskradius = 240;   % in pixels, AFTER binning operation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



% Create euler.txt file
eulerpath = fullfile(homepath,path,'euler.txt');
fid = fopen(eulerpath,'w');
fprintf(fid,'eu\n');
fprintf(fid,'%.0f\n',size(euler,1));
for ii = 1:size(euler,1)
    fprintf(fid,' %g %g %g\n',euler(ii,1),euler(ii,2),euler(ii,3));
end
fclose(fid);


% Create EMEBSD.nml file
nmlpath = fullfile(homepath,path,'EMEBSD.nml');
fid = fopen(nmlpath,'w');
fprintf(fid,'&EBSDdata\n');
fprintf(fid,' L = %.2f,\n',L);
fprintf(fid,' thetac = %.1f,\n',thetac);
fprintf(fid,' delta = %.1f,\n',delta);
fprintf(fid,' numsx = %g,\n',numsx);
fprintf(fid,' numsy = %g,\n',numsy);
fprintf(fid,' xpc = %.4f,\n',xpc);
fprintf(fid,' ypc = %.4f,\n',ypc);
fprintf(fid,' omega = %.1f,\n',omega);
fprintf(fid,' alphaBD = 0.0,\n');
fprintf(fid,' energymin = %.1f,\n',energymin);
fprintf(fid,' energymax = %.1f,\n',energymax);
fprintf(fid,' includebackground = ''y'',\n');
fprintf(fid,' anglefile = ''%seuler.txt'',\n',path);
fprintf(fid,' anglefiletype = ''orientations'',\n');
fprintf(fid,' eulerconvention = ''%s'',\n',eulerconvention);
fprintf(fid,' masterfile = ''%s'',\n',masterfile);
fprintf(fid,' energyfile = ''%s'',\n',masterfile);
fprintf(fid,' datafile = ''%sEBSDout.h5'',\n',path);
fprintf(fid,' bitdepth = ''8bit'',\n');
fprintf(fid,' beamcurrent = 150.0,\n');
fprintf(fid,' dwelltime = 100.0,\n');
fprintf(fid,' poisson = ''%s'',\n',poisson);
fprintf(fid,' binning = %g,\n',binning);
fprintf(fid,' applyDeformation = ''n'',\n');
fprintf(fid,' Ftensor = 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 1.D0,\n');
fprintf(fid,' scalingmode = ''%s'',\n',scalingmode);
fprintf(fid,' gammavalue = %g,\n',gammavalue);
fprintf(fid,' makedictionary = ''n'',\n');
fprintf(fid,' maskpattern = ''%s'',\n',maskpattern);
fprintf(fid,' maskradius = %g,\n',maskradius);
fprintf(fid,' hipassw = 0.05,\n');
fprintf(fid,' nregions = 10,\n');
fprintf(fid,' nthreads = %g,\n',nthreads);
fprintf(fid,' /\n');
fclose(fid);


% Run EMsoft EMEBSD
system(sprintf('EMEBSD %s',nmlpath));


% Extract images and save to file
h5path = fullfile(homepath,path,'EBSDout.h5');
data = h5read(h5path,'/EMData/EBSD/EBSDPatterns');
for ii = 1:size(data,3)
    img = data(:,:,ii)';    %% Use this line for EMsoft 4.3
    img = flipud(data(:,:,ii)');     % invert and flip array so that image is oriented properly (TSL convention)  %%Use this line for EMsoft 4.0, not for 4.3 
    s = strcat(outputname,'_',num2str(ii),imagetype);   % figure out file name
    outputpath = fullfile(homepath,path,s);
    imwrite(img,outputpath);
end


% delete files
[~,~] = system(sprintf('rm %s',eulerpath));
[~,~] = system(sprintf('rm %s',nmlpath));
[~,~] = system(sprintf('rm %s',h5path));

