function [euler, x, y, mad] = loadctf(inputpath, iconvert)
%LOADCTF
% Read in .ctf file outputted by EMFitOrientation
%%% Inputs:
% -inputpath: full path of .ctf file to read in
% -iconvert: convert euler angles to EMsoft convention? 1=yes, 0=no
%%% Outputs:
% -euler: Nx3 array containing phi1, PHI, phi2 (all in deg)
% -x: Nx1 array of x position
% -y: Nx1 array of y position
% -mad: Nx1 array of MAD
% 1/28/20 (Edward Pang, MIT)



% Read in data
fileID = fopen(inputpath);
C = textscan(fileID,'%s','Delimiter','\r');
fclose(fileID);

% figure out some things
datastartindex = 16;  % index of C where data begins
N = length(C{1}) - datastartindex + 1;  % number of map points
offset = datastartindex - 1;    % difference in array index between C and my stored data arrays below
    

% initialize some arrays
euler = zeros(N,3);
x = zeros(N,1);
y = zeros(N,1);
mad = zeros(N,1);


%%% Read in map data
% Loop through each data point
for ii=1:N
    % Extract data
    datarow = textscan(C{1}{ii+offset},'%f');   % data row delimited into cell array
    phi1 = datarow{1}(6);   % in deg
    PHI = datarow{1}(7);
    phi2 = datarow{1}(8);
    x(ii) = datarow{1}(2);
    y(ii) = datarow{1}(3);
    mad(ii) = datarow{1}(9);
    

    % Convert euler angles to EMsoft convention if specified
    %  phi1_oxford+90 = phi1_oemsoft (subtract 360 deg if necessary)
    %  PHI and phi2 are the same
    if iconvert==1
        phi1 = phi1+90;
        if phi1>=360
            phi1 = phi1-360;
        end
    end
    
    % store euler angles
    euler(ii,1) = phi1;
    euler(ii,2) = PHI;
    euler(ii,3) = phi2;
end


