function modelV = generateAllVoltages()
% Nate Faber (March 2014)
% Based loosely on code by Matt Johnson
% Program that reads the 3D coordinates of all the sections of the nucleus
% accumbens neurons. It then takes those 3D coordinates and calculates the
% maximum from DBS electrodes via a FEM. The voltages are then put into a
% hoc file which can be loaded by NEURON.
secPtsFileName = '../secPoints.txt';
femFileName = 'HumanDBS_axisymmetric.mph';
outputFileName = '../allVoltages.txt';

% Note - to get this matrix, you'll need to place four points in Rhino
% corresponding to [0,0,0], [0,0,3], c0, and c3.  Select c0 and c3 and
% create a block instance using the origin point of [0,0,0].  Then,
% rotate the c0/c3 points in the right view using [0,0,0] as the center
% of rotation.  Repeat for the front view as well.  Then translate the
% c0/c3 points until they match up with [0,0,0] and [0,0,3].  Click the
% 'Details' command button on the properties panel, and copy and then 
% paste the transformation matrix into the code below.
transmatrix = [0.88177443748849005, 0.10097511446662021, -0.46073622350571686, -1.5612511283791264e-17
        0.10097511446662021, 0.91375829790994301, 0.39350959233438881, 6.8087896432089678e-16
        0.46073622350571686, -0.39350959233438881, 0.79553273539843306, -6.8261368779687359e-16
        0, 0, 0, 1];

% Load the secPts file
secPts = importdata(secPtsFileName);
nPts = size(secPts,1);

% Transform the section coordinates into the correct place for the model
transPts = zeros(nPts,3);
for iPt = 1:nPts
    point = [secPts(iPt,:)*1e-3, 1]';
    newPoint = transmatrix * point;
    transPts(iPt,:) = newPoint(1:3)';
end

% Load the finite element model
model = mphload(femFileName);
model.hist.disable; % thought to free up memory issues with COMSOL v4.0

% Load the output file for writing
outFile = fopen(outputFileName, 'w');
    
%Convert 3D point into polar 2D
z = transPts(:,3)';
r = sqrt(sum(transPts(:,1:2) .^ 2,2))';
fprintf(1,'r: %f, %f\nz: %f, %f\n',min(r),max(r),min(z),max(z));

% Determine the COMSOL coordinates for cell compartments
modelV = mphinterp(model,'V2','coord',[r;z]);

%Write the value to an output file.
for iPt = 1:nPts
fprintf(outFile,'V_raw.x[%d] = %f\n',iPt-1,modelV(iPt));    
end

end