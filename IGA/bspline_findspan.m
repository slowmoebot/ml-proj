function i = bspline_findspan(n,p,xi,U)
% This routine searches for a given xi the corresponding knot span.
% n:        no. of control points
% p:        degree of the B-Spline/NURBS
% xi:       the xi we search the knot span for
% U(n+p+1): knot vector

%assert(xi >= 0 && xi<=1)

	% Check for xi = 1.
	if xi == 1
		i = n;
		return
	else
	% search for the knotspan i.
%  	i = find(U<=xi, 1, 'last' ); % backup function

	i=1;
	while(xi>=U(i))
		i=i+1;
	end
	i=i-1;
	%assert(i == find(U<=xi, 1, 'last' )) %check if correct
	end
end
