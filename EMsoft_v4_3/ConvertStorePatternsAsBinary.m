%CONVERTSTOREPATTERNSASBINARY
% Converts all images in a given folder to 8bit (if not already), saves 
% files with same name in a different location, also stores 8bit images 
% into a single .data file that can be read in by EMsoft
% Unlike EMsoft itself, this code can handle hexagonal and square grid data,
% automatically figures it out by reading in the .ang file.
% Reads in EDAX default file names in order of being acquired (start at top-
% left of map, row by row)
% ***File names must be sorted in order by MATLAB dir function, such as
% EDAX data that ends in the following: _00000, _00001, etc.
% 2/20/20 (Edward Pang, MIT) *adapted from some code from Marc De Graef

%%% Workaround for hexagonal grid data. Even numbered rows
% are shifted to the left by 0.5 x-spacings to make a square grid. Blank
% pattern is inserted wherever there is no data to fill a square grid. As a
% result, some of the parameters calculated by EMsoft may be incorrect. But
% the indexing is still ok. Programs in this EMsoft_utilities package knows
% this and unravels the data accordingly after indexing.


clear
tic

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
inputfolder = 'CuonCu/191031_647_Site3_lowmag/OIM Map 2/';  % Path (relative to EMData folder) to folder containing the input images
imagetype = '.png';     % This program will only look for images of this type (and will save it as the same type)
angfile = 'CuonCu/191031_647_Site3_lowmag/OIM Map 2/OIM Map 2.ang'; % Path (relative to EMData folder) to find the .ang file (Hough data)

%%% Outputs
% Path (relative to EMData folder) to folder where you want to save the 8bit images (you need to create it beforehand)
outputfolder = 'CuonCu/191031_647_Site3_lowmag/pattern_8bit/';

% Path (relative to EMData folder) of .data file to output
dataname = 'testdata2/scan_test.data';

% factor to bin patterns down by in .data file (=1 to keep current size)
%   will not bin down the individually outputted 8bit patterns
binning = 4;


%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



%%% error checking
% check if input and output folders are the same
if strcmp(inputfolder,outputfolder)
    error('Input and output folders are identical. Save output in a different locations.');
end

% check if datapath already exists
datapath = fullfile(homepath,dataname);     % full path
if exist(datapath,'file')
    error('.data file already exists. Rename dataname.');
end


%%% Figure out some things for .data file
% Read in grid info from .ang file
angpath = fullfile(homepath,angfile);
[~,~,grid] = loadang_header(angpath);

% figure out if hexagonal or square grid
if strcmp(grid{1},'HexGrid')
    issquare = 0;
else
    issquare = 1;
end

% get more grid info
ncols_odd = grid{4};
ncols_even = grid{5};
nrows = grid{6};


%%% find all images in folder
inputpath = fullfile(homepath,inputfolder);
outputpath = fullfile(homepath,outputfolder);
fileList = dir([inputpath '*' imagetype]);
N = length(fileList);   % how many files detected

% how many patterns should be present based on grid
if mod(nrows,2)==1  % if odd number of rows
    N_grid = (nrows-1)*(ncols_even+ncols_odd)/2 + ncols_odd;
else    % even number of rows
    N_grid = nrows*(ncols_even+ncols_odd)/2;
end
if N~=N_grid
    error('Number of patterns found (%g) does not match number of patterns expected (%g). Check folder and remove excess patterns.',N,N_grid);
end


%%% figure out some things from the first pattern
% see if patterns are already 8bit
inname = [inputpath fileList(1).name];   % Figure out full image path
im = imread(inname);
%     im = im(:,:,1);    % Extract first channel if multiple exist
if isa(im,'uint8')    % already uint8, directly save as .data file
    convertimg = 0;     % no need to convert
else
    convertimg = 1;     % needs to be converted to 8bit
end

% figure out image sizes
numsx = size(im,2);     % input
numsy = size(im,1);
imagewidth = numsx/binning;     % output
imageheight = numsy/binning;
  


%%% Main loop: convert patterns to 8bit (if necessary), output to .data file
fid = fopen(datapath,'w');       % Open .data file to write
nimg = 0;   % keep track of number of images outputted (including dummy patterns)

if convertimg==0    % already uint8, directly save as .data file
    % print warning
    warning('Patterns are already in 8bit format. Not outputting anything to: %s',outputfolder);
    
    % read in patterns and save to .data file
    for ii=1:N
        % Read in image
        inname = [inputpath fileList(ii).name];   % Figure out full image path
        im = imread(inname);
        im = im(:,:,1);    % Extract first channel if multiple exist
        im = single(im);    % Convert data to single precision
%         im = flipud(im);    % Flip image upside down (TSL has patterns upside down) [need this line for EMsoft v4.0, not for v4.3]

        % bin it down
        im_bin = zeros(imageheight,imagewidth,'uint8');   % initialize
        if binning>1
            for kk=0:binning:numsy-binning
                for ll=0:binning:numsx-binning
                    im_bin(kk/binning+1,ll/binning+1) = sum(sum(im(kk+1:kk+binning,ll+1:ll+binning)))/binning^2;
                end
            end
        else
            im_bin = im;
        end
        
        % Make a 1D vector of image
        p = zeros(imageheight*imagewidth,1);     % Initialize vector
        for jj = 1:imageheight
            for kk = 1:imagewidth
                p((jj-1)*imagewidth+kk) = im_bin(jj,kk);
            end
        end

        fwrite(fid,p,'float32');
        nimg = nimg + 1;
        
        
        % if at end of even row in hex grid, need to add a dummy pattern
        if issquare==0  % only an issue for hex grids
            if mod(ii,ncols_even+ncols_odd)==0
                im_bin = single(zeros(imageheight,imagewidth));     % Make dummy EBSD pattern (completely black) as missing data entry

                % Make a 1D vector of image
                p = zeros(imageheight*imagewidth,1);     % Initialize vector
                for jj = 1:imageheight
                    for kk = 1:imagewidth
                        p((jj-1)*imagewidth+kk) = im_bin(jj,kk);
                    end
                end

                fwrite(fid,p,'float32');
                nimg = nimg + 1;
            end
        end
        
        % print progress
        if mod(ii,1000)==0
            fprintf('  %g of %g completed.\n',ii,N);
        end
    end
    
    fprintf('%g images (%g dummy) saved in data file: %s\n\n',nimg,nimg-N,datapath);
else
    % Open each image, convert to 8bit, save to file and add to .data file
    for ii=1:N
        %%% Read in image
        inname = [inputpath fileList(ii).name];   % Figure out full image path
        im = imread(inname);
        im = im(:,:,1);    % Extract first channel if multiple exist

        
        %%% convert to uint8
        im = im2uint8(im);
        outname = [outputpath fileList(ii).name];
        imwrite(im,outname);
        
        
        %%% save to .data file
        im = single(im);    % Convert data to single precision
%         im = flipud(im);    % Flip image upside down (TSL has patterns upside down) [need this line for EMsoft v4.0, not for v4.3]

        % bin it down
        im_bin = zeros(imageheight,imagewidth,'uint8');   % initialize
        if binning>1
            for kk=0:binning:numsy-binning
                for ll=0:binning:numsx-binning
                    im_bin(kk/binning+1,ll/binning+1) = sum(sum(im(kk+1:kk+binning,ll+1:ll+binning)))/binning^2;
                end
            end
        else
            im_bin = im;
        end
        
        % Make a 1D vector of image
        p = zeros(imageheight*imagewidth,1);     % Initialize vector
        for jj = 1:imageheight
            for kk = 1:imagewidth
                p((jj-1)*imagewidth+kk) = im_bin(jj,kk);
            end
        end

        fwrite(fid,p,'float32');
        nimg = nimg + 1;
        
        
        % if at end of even row in hex grid, need to add a dummy pattern
        if issquare==0  % only an issue for hex grids
            if mod(ii,ncols_even+ncols_odd)==0
                im_bin = single(zeros(imageheight,imagewidth));     % Make dummy EBSD pattern (completely black) as missing data entry

                % Make a 1D vector of image
                p = zeros(imageheight*imagewidth,1);     % Initialize vector
                for jj = 1:imageheight
                    for kk = 1:imagewidth
                        p((jj-1)*imagewidth+kk) = im_bin(jj,kk);
                    end
                end

                fwrite(fid,p,'float32');
                nimg = nimg + 1;
            end
        end

        % print progress
        if mod(ii,1000)==0
            fprintf('  %g of %g completed.\n',ii,N);
        end
    end


    fprintf('\n%g patterns outputted to: %s\n',N,outputfolder);
    fprintf('%g images (%g dummy) saved in data file: %s\n\n',nimg,nimg-N,datapath);
end

fclose('all');



toc

