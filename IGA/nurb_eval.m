function S = nurb_eval(nurb, coeffs, dimS, points)

% Preallocate S
dim = numel(nurb.order);
% Sanity checks of the bspline structure
assert(dim == numel(nurb.number));
assert(size(points,1) == dim);
for i=1:dim
    assert(size(coeffs,i+1) == nurb.number(i))
end
assert(size(coeffs,1) >= dimS);

size_points = size(points,2);

p  = nurb.order(1) -1;

S = zeros(1,size_points);

% Perform two bspline_eval calls to determine the numerator and
% the denominator of Equation (4).
S = bspline_eval(nurb, coeffs, dimS, points);
S2 = bspline_eval(nurb, nurb.weights, 1, points);
% Divide coefficients by weights.


for j=1:size_points
    S(:,j)=S(:,j)./S2(1,j);
end

end
