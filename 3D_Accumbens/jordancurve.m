function [inout]=jordancurve(mesh,la)

% jordancurve() -- determine if a point is within a 3D mesh surface.  
%
% In this case, we will use the Jordan Curve Theorem, which states that a 
% point is within an object if the number of surface intersections from a 
% ray emanating from the point in an arbitrary direction is odd.  If the 
% number is even, the point is outside the object.
%
% Line: la + (lb-la)*t, where la=(xa,ya,za), lb=(xb,yb,zb)
% Plane: po + (p1-p0)*u+(p2-p0)*v, where pk=(xk,yk,zk), k=0,1,2
% Intersection point:
% [t]   [xa-xb x1-x0 x2-x0]-1[xa-x0]
% [u] = [ya-yb y1-y0 y2-y0]  [ya-y0]
% [v]   [za-zb z1-z0 z2-z0]  [za-z0]
% 
% If u,v are between [0,1] and (u+v)<=1, then the intersection point is in
% the plane inside the triangle spanned by pk [Moller-Trumbore, 1997].
%
% Inputs:
%  mesh - is a nx9 matrix with a row given by [x0,y0,z0,x1,y1,z1,x2,y2,z2]
%  la   - is the 3D coordinates of the point in question [xa,ya,za]
%
% Output:
%  inout - signal whether the point is inside (1) or outside (0) the mesh


% Declare variables
intersect_count=0;

% Define the ray segment projecting from (xa,ya,zy)
ray_dir = [0,0,50];      % arbitrary, change if needed
lb = la + ray_dir;

% Iterate through every triangle in the mesh
for i = 1:size(mesh,1),
    
    %Calculate the solution to the intersection point problem
    M_inv = inv([la(1)-lb(1), mesh(i,4)-mesh(i,1), mesh(i,7)-mesh(i,1); ...
        la(2)-lb(2), mesh(i,5)-mesh(i,2), mesh(i,8)-mesh(i,2); ...
        la(3)-lb(3), mesh(i,6)-mesh(i,3), mesh(i,9)-mesh(i,3)]);
    tuv = M_inv*[la(1)-mesh(i,1); la(2)-mesh(i,2); la(3)-mesh(i,3)];
    
    % Determine if the intersection point is in the plane inside the
    % triangle defined by pk; t must be [0,1] for the intersection point
    % to be on the ray segment
    if tuv(2)>=0 && tuv(2)<=1 && ...
            tuv(3)>=0 && tuv(3)<=1 && ...
            tuv(2)+tuv(3)<=1 && tuv(1)>=0,
        intersect_count = intersect_count + 1;
    end

end

% Determine if the point is within the 3D mesh surface
if mod(intersect_count,2)==0,
    inout = 0; % point outside
else
    inout = 1; % point inside
end