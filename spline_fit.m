function nurb = spline_fit(D,order,n_cp)

[n_points,n_dim]=size(D);

%n_cp = 4;
size_knot_vector = n_cp + order + 1;

crds = zeros(1,n_points-1);

for i = 1:n_points-1
    crds(i) = norm(D(i + 1, :) - D(i, :)); 
end

d = sum(crds);
xi_points = zeros(1, n_points);
for i = 1:n_points
    xi_points(i) = sum(crds(1:i-1)) / d;
end

xi = zeros(1, size_knot_vector);
xi(end - order:end) = 1;

xi(order+2:end-order)=(1:size_knot_vector-2*order-1)/(size_knot_vector-2*order-1);

dimS = 1;

bspline.number = n_cp;
bspline.knots{1} = xi;
bspline.order = order+1;
points = xi_points;

N = zeros(n_points,n_cp);
for i = 1:n_cp
    coeffs = zeros(1,n_cp);
    coeffs(i) = 1;
    N(:,i)=bspline_eval(bspline, coeffs, 1, points);
end

P = N\D;

nurb.knots{1} = xi;
nurb.weights = ones(1, n_cp);
nurb.coords = P';
nurb.number = n_cp;
nurb.order = order+1;

end