function draw_nurb_surf(nurb,points, fhandle)

% nurb:      NURBS data structure
% points(2): the number of subintervals of [0,1] in one direction,
%            used to draw the NURBS as a set of quadrilaterals
% fhandle:   Function that will be plotted on top of the surface as a
%            contour plot, e.g. @(x,y) x^2

assert(numel(points) == 2);
assert(isa(fhandle, 'function_handle'));

deltaX = 1/(points(1)-1);
deltaY = 1/(points(2)-1);

[X, Y] = meshgrid(0:deltaX:1.0, 0:deltaY:1.0);

points_2d = zeros(2, points(1)*points(2));
points_2d(1,:) = reshape(X,1,[]);
points_2d(2,:) = reshape(Y,1,[]);

S = nurb_eval(nurb,nurb.coords,2,points_2d);

for i=1:points(2)
    for j=1:points(1)
        Z(i,j) = fhandle(X(i,j),Y(i,j));
    end
end

Xsurf = reshape(S(1,:), points(2), points(1));
Ysurf = reshape(S(2,:), points(2), points(1));


s = surf(Xsurf, Ysurf, Z);
s.FaceColor = 'interp';
c = colorbar;
colormap('jet');
% proper aspectratio
daspect([ 2 2 1 ]);
% topview
view(0,90);


drawnow;

end
