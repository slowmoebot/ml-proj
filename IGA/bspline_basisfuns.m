function N = bspline_basisfuns(i,xi,p,U)
% This routine determines for a given xi and its corresponding knot span i
% the non-vanishing basis functions N(1:p+1), which are evaluated at xi.
% i:        to xi corresponding knot span
% xi:       the xi we search the knot span for
% p:        degree of the B-Spline/NURBS
% U(n+p+1): knot vector

% Preallocate the N array
	N = zeros(p+1,1);

	% First iteration of the outer loop
	
	N(1) = 1;
		
	
	for j=2:p+1 % iterates over columns 
		
		% For k=1 there is no dependence on N(k-1) of the previous run.
		saved = 0.0;
		
		for k=1:j-1 % iterates over rows
			
			% Compute $N_{i-j+k,j-1}$ according to the Cox-deBoor formula
			
			tmp = saved +  (U(i+k) - xi) /(U(i+k) - U(i-j+k+1)) * N(k);
			

			% Precompute the first summand of $N_{i-j+k+1,j-1}$ for the next
			% iteration of the inner loop.
			
			saved = (xi - U(i-j+k+1) ) / (U(i+k) - U(i-j+k+1)) * N(k);
			
			% Assign to basis function vector
			N(k) = tmp;
			
		end
		
		% Unroll the last iteration of the inner loop 
		N(j) = saved;
		
	end
	%sum(N)
	%check if partition of unity holds
	assert(abs(1.0 - sum(N)) <= eps(1))

end
