function internalPoints = generateInternalPoints(mesh,nPoints)
%Uses the Jordan Curve Algorithm to generate a list of points that
%fall within a given mesh

maxmat = max(mesh,[],1);
minmat = min(mesh,[],1);

xmin = min(minmat([1,4,7]));
xmax = max(maxmat([1,4,7]));
ymin = min(minmat([2,5,8]));
ymax = max(maxmat([2,5,8]));
zmin = min(minmat([3,6,9]));
zmax = max(maxmat([3,6,9]));

iPoint = 1;
internalPoints = zeros(nPoints,3);

drawMesh(mesh);
while iPoint <= nPoints
    xpos = random('unif',xmin,xmax);
    ypos = random('unif',ymin,ymax);
    zpos = random('unif',zmin,zmax);
    if jordancurve(mesh,[xpos,ypos,zpos])
        internalPoints(iPoint,:) = [xpos, ypos, zpos];
        scatter3(xpos,ypos,zpos,'r');
        drawnow;
        if mod(iPoint,10) == 0
            disp(iPoint);
        end
        iPoint = iPoint + 1;
    end
end