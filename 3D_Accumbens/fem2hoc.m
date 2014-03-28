% fem2hoc
% Nate Faber (March 2014)
% Based loosely on code by Matt Johnson
% Program that reads the 3D coordinates of all the sections of the nucleus
% accumbens neurons. It then takes those 3D coordinates and calculates the
% maximum from DBS electrodes via a FEM. The voltages are then put into a
% hoc file which can be loaded by NEURON.

secPtsFile = '../secPts.txt';

% Note - to get this matrix, you'll need to place four points in Rhino
% corresponding to [0,0,0], [0,0,3], c0, and c3.  Select c0 and c3 and
% create a block instance using the origin point of [0,0,0].  Then,
% rotate the c0/c3 points in the right view using [0,0,0] as the center
% of rotation.  Repeat for the front view as well.  Then translate the
% c0/c3 points until they match up with [0,0,0] and [0,0,3].  Click the
% 'Details' command button on the properties panel, and copy and then 
% paste the transformation matrix into the code below.
transmatrix = [0.88177443748849005, 0.10097511446662021, -0.46073622350571686, -11.065361701000583;
    0.10097511446662021, 0.91375829790994301, 0.39350959233438881, -3.0420521151039122;
    0.46073622350571686, -0.39350959233438881, 0.79553273539843306, 3.3673465819877424;
    0, 0, 0, 1];

% Transform the cell population coordinate space to COMSOL
numcells = size(pop,2);
h = waitbar(0,'Transforming cell population...');
for i = 1:numcells,
    
    
    

    waitbar(i/size(pop,2),h,['Cell: ',num2str(i)]);
    
    % Determine the COMSOL coordinates for cell compartments
    for j = 1:size(pop(i).pts,1),
        popok(i).pts(j,:) = pop(i).pts(j,:) * transmatrix(1:3,1:3)' + transmatrix(:,4)';
    end
    
end
close(h);

clearvars -except popok model filename numcells

for i=1:numcells,
    popok(i).vpts = zeros(size(popok(i).pts,1),size(popok(i).pts,2));
end

% Determine the voltage at each compartment
model.hist.disable; % thought to free up memory issues with COMSOL v4.0
h=waitbar(0,'Calculating voltages...');
for k = 1:size(popok,2),        % number of cells
    waitbar(k/size(popok,2),h,['Cell: ',num2str(k)]);
    popok(k).vpts = mphinterp(model,'V','Coord',[popok(k).pts(:,1),popok(k).pts(:,2),popok(k).pts(:,3)]');
end
close(h);

% Create a text file containing .txt code of the segment voltages
for i = 1:numcells,   % print compartments
    ctr = 0;
    fid = fopen([filename,num2str(i),'.txt'],'w');
    for j = 1:length(popok(i).vpts),  
        fprintf(fid,'%s%f\n',['V_raw.x[',num2str(ctr),']='],popok(i).vpts(j));
        ctr = ctr + 1;
    end
    fclose(fid);
end

toc