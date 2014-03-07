function drawMesh(meshMat,linespec)
if nargin < 2
    linespec = 'b';
end
hold on;
for i = 1:size(meshMat,1)
    xPoints = meshMat(i,[1, 4, 7, 1]);
    yPoints = meshMat(i,[2, 5, 8, 2]);
    zPoints = meshMat(i,[3, 6, 9, 3]);
    plot3(xPoints,yPoints,zPoints,linespec);
end