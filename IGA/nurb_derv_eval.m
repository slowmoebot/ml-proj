function dS = nurb_derv_eval(nurb, coeffs, dimS, points)

dim = numel(nurb.order);
assert(dim == numel(nurb.number));
assert(size(points,1) == dim);
for i=1:dim
	assert(size(coeffs,i+1) == nurb.number(i))
end
assert(size(coeffs,1) >= dimS);

size_points = size(points,2);

% Preallocate dS
dS = zeros(dimS,dim,size_points);

% Compute A,W and the corresponding jacobians dA,dW
size_points = size(points,2);
% 
% p = nurb.order(1) -1;
% q = nurb.order(2) -1;
% 
% S = zeros(1,size_points);

%coeffs = ones(1,nurb.number(1),nurb.number(2)) .* nurb.weights;

% Perform two bspline_eval calls to determine the numerator and
% the denominator of Equation (4).

A = bspline_eval(nurb, coeffs, dimS, points);			% 2 x sizepoints
W = bspline_eval(nurb, nurb.weights, 1, points);		% 1 x sizepoints

dA = bspline_derv_eval(nurb, coeffs, dimS, points);		% 2 x 2 x sizepoints
dW = bspline_derv_eval(nurb, nurb.weights, 1, points);  % 2 x 1 x sizepoints

for k=1:size_points
	% Compute dS(:,:,k)
	
	dS(:,:,k) = 1 / W(:,k) * (dA(:,:,k) - A(:,k) * dW(:,:,k) / W(:,k));
end
end



%dS(:,:,k)=(dA(:,:,k) - A(:,k)./W(:,k).*W(:,k))./W(:,k);