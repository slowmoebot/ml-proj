function S = bspline_eval(bspline, coeffs, dimS, points)

% bspline:              holds the geometrical B-Spline structure
% coeffs(dimS,n1[,n2]): the coefficients we should use for evaluation, which are NOT
%                       necessarily the control points
% dimS:                 dimension of the space the curve/surface/DOF is embedded in,
%                       which implies that the first index of coeffs should have at
%                       least that many entries.
% points(dim,size_points): the evaluation points with points(1,i) = $\xi_i$ and points(2,i) = $\eta_i$

dim = numel(bspline.order);
% Sanity checks of the bspline structure
assert(dim == numel(bspline.number));
assert(size(points,1) == dim);
for i=1:dim
	assert(size(coeffs,i+1) == bspline.number(i))
end
assert(size(coeffs,1) >= dimS);

size_points = size(points,2);

%Preallocate S
%YOUR_CODE

S = zeros(dimS,size_points);



if (dim == 1)
	
	% degree of bspline
	p = bspline.order-1;
	% number of control points
	n = bspline.number;
	for j=1:size_points
		% For each evaluation point points(:,j) find the corresponding knotspan i
		% using bspline_findspan and then determine the $p+1$ non-zero
		% basis functions N with bspline_basisfuns. Use this information to
		% form the linear combination.
		xi = points(j);
		U = bspline.knots{1};
		i = bspline_findspan(n,p,xi,U);
		N = bspline_basisfuns(i,xi,p,U);
		
		for k=1:p+1
			for d=1:dimS
				S(d,j) = S(d,j) + coeffs(d,i-p+k-1) * N(k);
			end
		end
		
	end

elseif (dim == 2)
	
	p  = bspline.order(1) -1;
	q  = bspline.order(2) -1;
	
	n1 = bspline.number(1);
	n2 = bspline.number(2);
	
	U1 = bspline.knots{1};
	U2 = bspline.knots{2};
	
	
	for j=1:size_points
		% For each evaluation point points(:,j) find the corresponding knotspan i
		% using bspline_findspan and then determine the $p+1$ non-zero
		% basis functions N with bspline_basisfuns. Use this information to
		% form the linear combination.
		xi = points(1,j);
		eta = points(2,j);
		
		
		i1 = bspline_findspan(n1,p,xi,U1);
		i2 = bspline_findspan(n2,q,eta,U2);
		
		N1 = bspline_basisfuns(i1,xi,p,U1);
		N2 = bspline_basisfuns(i2,eta,q,U2);
		
		for k1=1:p+1
			for k2=1:q+1
				for d=1:dimS
					
					S(d,j) = S(d,j) + coeffs(d,i1-p+k1-1,i2-q+k2-1) * N1(k1) * N2(k2);

			end
		end
	end
	
end

end



