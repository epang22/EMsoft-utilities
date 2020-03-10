function [PC, phaseinfo, grid] = loadang_header(inputpath)
%LOADANG
% Read in header data from .ang file
%%% Inputs:
% -inputpath: full path of .ang file to read in
%%% Outputs:
% -PC: [X*, Y*, Z*] in fraction detector width
% -phaseinfo: 1x4 cell array. each row is phase. col1=phase id, col2=materialname, col3=symmetry id, col4: lattice param [a b c alpha beta gamma] (A, deg)
% -grid: 1x6 cell array [gridtype xstep ystep ncols_odd ncols_even nrows]
% 1/16/20 (Edward Pang, MIT)



% Read in data
fileID = fopen(inputpath);
C = textscan(fileID,'%s','Delimiter','\r');
fclose(fileID);


% Loop through each line to find start of map data
for ii=1:length(C{1})
    if ~isempty(regexp(C{1}{ii},'SCANID','once'))
        datastartindex = ii+2;  % index of C where data begins
        break
    end
end


% Read in PC info
for ii=1:length(C{1})
    if ~isempty(regexp(C{1}{ii},'x-star', 'once'))
        dataxstar = textscan(C{1}{ii},'%s %s %f');
        xstar = dataxstar{3};
    end
    if ~isempty(regexp(C{1}{ii},'y-star', 'once'))
        dataystar = textscan(C{1}{ii},'%s %s %f');
        ystar = dataystar{3};
    end
    if ~isempty(regexp(C{1}{ii},'z-star', 'once'))
        datazstar = textscan(C{1}{ii},'%s %s %f');
        zstar = datazstar{3};
        break
    end
end
PC = [xstar ystar zstar];   % put in single vector for output


% Figure out phases, their symmetry, and their lattice params
phaseinfo = {};  % each row is phase, {1}=phase id, {2}: material name, {3}=symmetry id, {4}: lattice param [a b c alpha beta gamma] (A, deg)
counter = 0;    % number of phases
for ii=1:length(C{1})
    np = regexp(C{1}{ii},'Phase','once');    % if not empty, equals starting index of 'Phase' in string
    nn = regexp(C{1}{ii},'MaterialName','once');   % if not empty, equals starting index of 'MaterialName' in string
    ns = regexp(C{1}{ii},'Symmetry','once');    % if not empty, equals starting index of 'Symmetry' in string
    nlp = regexp(C{1}{ii},'LatticeConstants','once');    % if not empty, equals starting index of 'Symmetry' in string
    
    % get phase id
    if ~isempty(np)
        counter = counter+1;
        wholestring = C{1}{ii};
        phaseinfo{counter,1} = str2double(wholestring(np+6:end));  % store phase ID
    end

    % get material name
    if ~isempty(nn)
        wholestring = C{1}{ii};
        tabindex = regexp(wholestring,'\t','once');
        phaseinfo{counter,2} = wholestring(tabindex+1:end);  % store material name
    end
    
    % get symmetry id
    if ~isempty(ns)
        wholestring = C{1}{ii};     % whole line as single string
        wholestring2 = textscan(wholestring,'%s %s %f');    % space delimited into cell array
        phaseinfo{counter,3} = wholestring2{3};
    end

    % get lattice parameters 
    if ~isempty(nlp)
        wholestring = C{1}{ii};     % whole line as single string
        wholestring2 = textscan(wholestring,'%s %s %f %f %f %f %f %f');    % space delimited into cell array
        phaseinfo{counter,4} = [wholestring2{3} wholestring2{4} wholestring2{5} wholestring2{6} wholestring2{7} wholestring2{8}];
    end

    % stop once past header region
    if ii>datastartindex
        break
    end
end



% Read in grid info
for ii=1:length(C{1})
    if ~isempty(regexp(C{1}{ii},'GRID', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %s');
        grid{1} = gridtemp{3}{1};
    end
    if ~isempty(regexp(C{1}{ii},'XSTEP', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{2} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'YSTEP', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{3} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'NCOLS_ODD', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{4} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'NCOLS_EVEN', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{5} = gridtemp{3};
    end
    if ~isempty(regexp(C{1}{ii},'NROWS', 'once'))
        gridtemp = textscan(C{1}{ii},'%s %s %f');
        grid{6} = gridtemp{3};
        break
    end
end


