%COMPAREPATTERNS
% Load 2 EBSD patterns, compute dot product, plot intensity difference
% and a window that flashes between the 2 patterns
% Fill in INPUT PARAMETERS section with desired parameters
% 1/13/20 (Edward Pang, MIT)

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path1 = 'testdata2/Scan2_1.png';  % path within EMdatapathname to pattern1
path2 = 'testdata2/Scan2_1_sim.png';  % path within EMdatapathname to pattern2

timemax = 10;   % how many seconds to flash between the two patterns


%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



% load experimental patterns
fullpath1 = fullfile(homepath,path1);
fullpath2 = fullfile(homepath,path2);

% remove extra channels if necessary
img1 = imread(fullpath1);
if size(img1,3)>1
    img1 = img1(:,:,1);
end
img2 = imread(fullpath2);
if size(img2,3)>1
    img2 = img2(:,:,1);
end

% check that images are the same size
if size(img1,1)~=size(img2,1) || size(img1,2)~=size(img2,2)
    error('Images are not the same size!');
end

% dot product computations
img1norm = double(img1);
img1norm = img1norm./sqrt(sum(sum(img1norm.^2)));
img2norm = double(img2);
img2norm = img2norm./sqrt(sum(sum(img2norm.^2)));
dp = sum(sum(img1norm.*img2norm));

% convert images to uint8 if not already
if ~isa(img1,'uint8')
    img1 = im2uint8(img1);
end
if ~isa(img2,'uint8')
    img2 = im2uint8(img2);
end

% difference plot
imgdiff = abs(img1 - img2);
diffmax = max(max(imgdiff));
imshow(imgdiff,'Colormap',hot(255));
set(gcf,'Position',[100 100 size(imgdiff,1) size(imgdiff,2)]);
set(gca,'Position',[0 0 1 1],'units','normalized');
colorbar

% print info
fprintf('Dot product = %.6f\n',dp);
fprintf('Max difference = %g\n',diffmax);

% flash between two images
figure(2);
for ii=1:timemax
    imshow(img1); pause(0.5);
    imshow(img2); pause(0.5);
end
