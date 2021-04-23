%FINDMAPPOINTS
% Load in data from .ang file, plot map, click points on map to 
% obtain coordinates and other info
% EDAX/TSL coordinate system
% Fill in INPUT PARAMETERS section with desired parameters
% Original: 2/22/20 (Edward Pang, MIT)
% Change log:
% -4/23/21 ELP: fix bug was outputting 1 row off, change output to 
% 0-indexing, add phase to output

clear

%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path (relative to EMData folder) to find the .ang file
angfile = 'testdata2/Scan1.ang';

plotoption = 5;     % what variable to color map with (1=phi1, 2=PHI, 3=phi2, 4=CI, 5=IQ, 6=fit)
markersize = 3;     % use smaller value if you have more data points (3 is a reasonable value for 20000 map points)


%%% Parameters you don't need to change often %%%
% Paths for this computer
homepath = '/home/eddie/EMsoftfiles/EMData/';  % path to EMdatapathname (don't need to touch after initial setup of EMsoft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW THIS LINE (unless you are doing code development) %



% Read in .ang file
fullpath = fullfile(homepath,angfile);
[euler, x, y, IQ, CI, phase, fit, ~, grid, ~] = loadang(fullpath, 1);

N = length(x);  % number of map points
xstep = grid{2};
ystep = grid{3};
ncolsodd = grid{4};
ncolseven = grid{5};
nrows = grid{6};
phi1 = euler(:,1);
PHI = euler(:,2);
phi2 = euler(:,3);



% pick one of these quantities for plotting
if plotoption==1
    data = phi1;
elseif plotoption==2
    data = PHI;
elseif plotoption==3
	data = phi2;
elseif plotoption==4
    data = CI;
elseif plotoption==5
    data = IQ;
elseif plotoption==6
    data = fit;
else
    error('Invalid plotoption input.');
end

% vector of x and y values
xvectorodd = 0:xstep:xstep*(ncolsodd-1);
xvectoreven = xstep/2:xstep:(xstep*(ncolseven-1)+xstep/2);
yvector = 0:ystep:(nrows-1)*ystep;

% plot
Ncolors = 256;  % number of colors to discretize into (higher=prettier image, slower plotting)
colors = parula(Ncolors);   % rgb of colors in each row
datamin = min(data);    % min of data for color scaling
datarange = max(data)-min(data);    % range of data for color scaling


% create figure window
figure('Position',[100 100 800 round((max(y)/max(x))*800)]);
h = axes;

% loop through each color, figure out corresponding data, add to plot
for ii=1:Ncolors
    % get bounding values of data for this color
    p_low = (ii-1)/Ncolors;     % lower percentile of data to plot this iter
    p_high = ii/Ncolors;        % higher percentile of data to plot this iter
    data_low = p_low*datarange + datamin;  % min value of data to plot this iter
    data_high = p_high*datarange + datamin;  % max value of data to plot this iter
    
    % extract data for this color
    index = (data>=data_low) & (data<data_high);     % indices of points to plot
    
    % add to plot
    plot(x(index),y(index),'o','MarkerEdgeColor',colors(ii,:),'MarkerFaceColor',colors(ii,:),'MarkerSize',markersize,...
        'ButtonDownFcn',{@ImageClickCallback,xvectorodd,xvectoreven,yvector,ncolsodd,ncolseven,phi1,PHI,phi2,CI,IQ,fit,phase}); hold on;
end

% add max point
[~,imax] = max(data);
plot(x(imax),y(imax),'o','MarkerEdgeColor',colors(end,:),'MarkerFaceColor',colors(end,:),'MarkerSize',markersize,...
    'ButtonDownFcn',{@ImageClickCallback,xvectorodd,xvectoreven,yvector,ncolsodd,ncolseven,phi1,PHI,phi2,CI,IQ,fit,phase});

% tweak plot
set(gca,'Ydir','reverse');  % OIM map convention y pointing down
set(h, 'Color', 'None', 'Xtick', [], 'Ytick', []);
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);


% print column labels
fprintf('Index numbers begin with 0 for the first pattern.\n');  
fprintf(' index:        x,        y,       phi1,        PHI,       phi2,       CI,         IQ,      fit,    phase\n');  



% define function for click callback
function ImageClickCallback(objectHandle,~,xvectorodd,xvectoreven,yvector,ncolsodd,ncolseven,phi1,PHI,phi2,CI,IQ,fit,phase)
    % extract x and y coordinates where you clicked    
    axesHandle = get(objectHandle,'Parent');
    coordinates = get(axesHandle,'CurrentPoint'); 
    xclick = coordinates(1,1);
    yclick = coordinates(1,2);
    
    % figure out nearest point to where you clicked
    [~,iy] = min(abs(yvector-yclick));
    y = yvector(iy);    % nearest y value
    
    if mod(iy,2)==0   % even row
        [~,ix] = min(abs(xvectoreven-xclick));
        x = xvectoreven(ix);    % nearest x value
        index = (iy/2)*ncolsodd + (iy/2-1)*ncolseven + ix;  % linear index of this point (1-index)
    else    % odd row
        [~,ix] = min(abs(xvectorodd-xclick));
        x = xvectorodd(ix);    % nearest x value
        index = ((iy-1)/2)*ncolsodd + ((iy-1)/2)*ncolseven + ix;  % linear index of this point (1-index)
    end
    
    hold on; plot(xclick,yclick,'ok');    % plot marker at nearest grid point to where you clicked
    
    % print relevant info to screen
    fprintf('%6.0f: %8.4f, %8.4f, %10.4f, %10.4f, %10.4f, %8.4f, %10.4f, %8.4f, %8.0f\n',index-1,x,y,phi1(index),PHI(index),phi2(index),CI(index),IQ(index),fit(index),phase(index));
end

