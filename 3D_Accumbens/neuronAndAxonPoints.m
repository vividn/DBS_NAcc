function neuronAndAxonPoints()
% Make a list of points for neurons in the nucleus accumbens and a
% corresponding list of points in the ventral pallidum. Points are grouped
% by ascending x-value so that medial neurons will have medial axons.
close all;
clear all;
warning('off','MATLAB:singularMatrix')

nNeurons = 100;
nDends = 4;
somaLength = 0.015; %in mm
dendLengthArray = [20, 24.23, 395.2]./1000; % from Wolf et al 2005
% Final number of points in each axon contour
ptsPerAxon = 100;

Acc = importdata('AccMesh.raw');
Acc = Acc.data;

VP = importdata('VPMesh.raw');
VP = VP.data;

figure(1)
%Plot DBS electrode line
DBSx = [7.51287;12.3506];
DBSy = [8.22211;4.09026];
DBSz = [-6.57997;1.77312];

plot3(DBSx,DBSy,DBSz,'m','LineWidth',4);

accPoints = generateInternalPoints(Acc,nNeurons);
accPoints = sortrows(accPoints,1);

vpPoints = generateInternalPoints(VP,nNeurons);
vpPoints = sortrows(vpPoints,1);

close all;
figure(1);
hold on;
scatter3(accPoints(:,1),accPoints(:,2),accPoints(:,3),'g');
scatter3(vpPoints(:,1),vpPoints(:,2),vpPoints(:,3),'r');
plot3(DBSx,DBSy,DBSz,'m','LineWidth',4);

% The accPoints are where the neuron and the axon connect
% The vpPoints are where the axon ends.

% Generate the axon contours to procedurally generate different axon
% sections:
allAxonPts = NaN(ptsPerAxon,3,nNeurons);
allSomaPts = NaN(2,3,nNeurons);
allDends = cell(nDends,nNeurons);
for iNeuron = 1:nNeurons
    axonPts = directedRandomWalk(accPoints(iNeuron,:),vpPoints(iNeuron,:),.4,.8);
    nAxonPts = size(axonPts,1);
    
    % Smooth axon using spline
    step = (nAxonPts-1)/(ptsPerAxon-1);
    F = spline(1:nAxonPts,axonPts');
    smoothAxon = ppval(F,1:step:nAxonPts)';
    allAxonPts(:,:,iNeuron) = smoothAxon;
    
    % Make soma in same direction as first axon segment
    somaDir = unitvec(smoothAxon(1,:) - smoothAxon(2,:));
    somaPts = [smoothAxon(1,:) + somaDir*somaLength; smoothAxon(1,:)];
    allSomaPts(:,:,iNeuron) = somaPts;
    
    % Create dendrites
    for iDend = 1:nDends
        if (iDend+1)/nDends <= 0.5
            somaSide = 1;
            somaVec = somaDir*1;
        else
            somaSide = 2;
            somaVec = somaDir*-1;
        end
        allDends{iDend,iNeuron} = branchingDend(somaPts(somaSide,:),somaVec,dendLengthArray);
    end
    
    plot3(smoothAxon(:,1),smoothAxon(:,2),smoothAxon(:,3),'b'); 
end


axon_pts2hoc(allAxonPts,allSomaPts,allDends);
end

function DendStruct = branchingDend(startPoint,inputVector,lenArray)
% Creates a branching dendrite recursively
% Creates two dendrite pieces going in random direction, but within a
% certain angle of the input segment
% (the projection of the new vectors on the inputVector is greater than
% cos(pi/4)).
% The new dendrite pieces have the length of the first entry in lenArray
% The process is then repeated recursively gradually shortening lenArray
% until it is empty.
%
% Output: DendStruct
%           -> startPoint 
%           -> endPoint
%           -> child1 (if it exists)
%           -> child2 (if it exists)

DendStruct.startPoint = startPoint;
proj = -1;
while proj < cos(pi/4)
    % choose a random direction and see if forward movement occurs
    newVector = unitvec(rand(1,3)-0.5);
    proj = dot(newVector,unitvec(inputVector));
end
endPoint = startPoint + newVector * lenArray(1);
plot3([startPoint(1);endPoint(1)],[startPoint(2);endPoint(2)],[startPoint(3);endPoint(3)],'g')
DendStruct.endPoint = endPoint;
if length(lenArray) > 1
    DendStruct.child1 = branchingDend(endPoint,newVector,lenArray(2:end));
    DendStruct.child2 = branchingDend(endPoint,newVector,lenArray(2:end));
end
end



    