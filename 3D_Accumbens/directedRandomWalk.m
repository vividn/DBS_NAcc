function pathVerts = directedRandomWalk(startPoint,targetPoint,stepSize,randomness)
% Function creates a random path from the startPoint to the targetPoint.
% The program uses a biased random walk, which adds a random direction
% vector with a vector always pointing towards the target. The vector
% pointing to the target is of magnitude stepSize, and the random vector is
% of size randomness*stepsize.
%

distToTarget = sqrt(sum((targetPoint - startPoint).^2));
totalDist = distToTarget;

%scatter3(startPoint(1), startPoint(2), startPoint(3),90,'filled','g');
%scatter3(targetPoint(1), targetPoint(2), targetPoint(3),90,'filled','r');
if nargin < 4
    randomness = 1;
    if nargin < 3
        stepSize = distToTarget / 100;
    end
end

pathVerts = zeros(round(10*distToTarget/stepSize),3);
pathVerts(1,:) = startPoint;
lastPoint = startPoint;
nVerts = 1;

while distToTarget > stepSize*randomness + stepSize
    % randomness decreases while approaching target
    rndFactor = randomness * distToTarget / totalDist;
    randvec = rndFactor * stepSize * unitvec(rand(1,3)-0.5);
    biasvec = stepSize * unitvec(targetPoint - lastPoint);
    lastPoint = lastPoint + randvec + biasvec;
    nVerts = nVerts+1;
    pathVerts(nVerts,:) = lastPoint;
    distToTarget = sqrt(sum((targetPoint - lastPoint).^2));
    %scatter3(lastPoint(1),lastPoint(2),lastPoint(3),30,'filled','b');
    %drawnow;
end
nVerts = nVerts+1;
pathVerts(nVerts,:) = targetPoint;
pathVerts = pathVerts(1:nVerts,:); %truncate extra slots
%plot3(pathVerts(:,1),pathVerts(:,2),pathVerts(:,3),'--b');

end