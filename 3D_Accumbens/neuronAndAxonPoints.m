% Make a list of points for neurons in the nucleus accumbens and a
% corresponding list of points in the ventral pallidum. Points are grouped
% by ascending x-value so that medial neurons will have medial axons.
close all;
clear all;
warning('off','MATLAB:singularMatrix')

nNeurons = 20;
somaLength = 0.02; %in mm
% Final number of points in each axon contour
ptsPerAxon = 1000;

Acc = importdata('AccMesh.raw');
Acc = Acc.data;

VP = importdata('VPMesh.raw');
VP = VP.data;

figure(1)
accPoints = generateInternalPoints(Acc,nNeurons);
accPoints = sortrows(accPoints,1);
clf;

vpPoints = generateInternalPoints(VP,nNeurons);
vpPoints = sortrows(vpPoints,1);

close all;
figure(1);
hold on;
scatter3(accPoints(:,1),accPoints(:,2),accPoints(:,3),'g');
scatter3(vpPoints(:,1),vpPoints(:,2),vpPoints(:,3),'r');

% The accPoints are where the neuron and the axon connect
% The vpPoints are where the axon ends.

% Generate the axon contours to procedurally generate different axon
% sections:
allAxonPts = NaN(ptsPerAxon,3,nNeurons);
allSomaPts = NaN(2,3,nNeurons);
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
    
    plot3(smoothAxon(:,1),smoothAxon(:,2),smoothAxon(:,3),'b'); 
end
axon_pts2hoc(allAxonPts,allSomaPts);
    