%GETPATTERNCENTER_EDAXANG
% Get pattern center in EMsoft coordinates from EDAX .ang file
% Fill in INPUT PARAMETERS section with desired parameters
% Original: 1/16/20 (Edward Pang, MIT)
% Change log:
% -4/24/21 ELP: change xpc sign to match new coordinate system

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path (relative to EMData folder) to find the .ang file (Hough data)
angfile = 'testdata2/Scan1.ang';


%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)

% Detector parameters
delta = 7.4*480/416;   % CCD pixel size on the scintillator surface (microns)
numsx = 416;    % number of CCD pixels along x and y
numsy = 416;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



% read in ang file
fullpath = fullfile(homepath,angfile);
PC = loadang_header(fullpath);
xstar = PC(1);
ystar = PC(2);
zstar = PC(3);


% convert to EMsoft coordinates
L = zstar*numsx*delta;
xpc = -1*numsx*(xstar-0.5);
ypc = numsy*(ystar-0.5);


% print
fprintf('Pattern center (EMsoft coordinates):\n');
fprintf('  xpc = %.4f px\n',xpc);
fprintf('  ypc = %.4f px\n',ypc);
fprintf('  L = %.2f um\n',L);

