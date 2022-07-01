function [mat rhs] = assemble_matrix_2d(use_bspline, ndof, nurb, ...
    element_matrix, rhs_function)

%% Preamble (exercise 4)

% preamble and quadrature part of the code
% don't worry about this for the moment
assert(isa(element_matrix, 'function_handle'));
assert(isa(use_bspline, 'logical'));

% number of shape functions per element
nshp = nurb.order(1) * nurb.order(2);
alldof = ndof*nurb.number(1)*nurb.number(2);
mat = spalloc(alldof,alldof,alldof*nshp);
rhs = zeros(alldof,1);

%remove double entries in the knot arrays
xknots = unique(nurb.knots{1});
yknots = unique(nurb.knots{2});

[quad_weight_1d, quad_points_1d] = quad_rule();
nquad_1d = size(quad_weight_1d, 2);

nquad = nquad_1d * nquad_1d;
quad_weight = zeros(1,nquad);
quad_points = zeros(2,nquad);

nurb_area = 0.0;

% Loop over all elements, which are spaced from knot to knot
for ie = 1:size(xknots,2)-1
	for je = 1:size(yknots,2)-1
		
		%% Exercise 4
		
        % build quad rules
        xi1 = xknots(ie);
        eta1 = yknots(je);
        deltaXi = xknots(ie+1) - xi1;
        deltaEta = yknots(je+1) - eta1;
        % area of the element in the reference domain
        area = deltaXi * deltaEta;
        % Compute the weights according to equation (4).
        % It is easier if you use the function 'reshape', which makes
        % it possible to avoid the integer division and the modulo operation.
		
		quad_weight = reshape(quad_weight_1d' * quad_weight_1d,1,nquad) / 4 * area;
		
        % Compute the quadrature points according to formula (5).
        % If possible use the functions 'repmat' and 'reshape'.
		
        quad_points(1,:) = reshape(repmat(quad_points_1d  * deltaXi / 2,nquad_1d,1),1,nquad);
		quad_points(1,:) = quad_points(1,:) + repmat(xi1 + deltaXi / 2,1,nquad);
		
        quad_points(2,:) = repmat(quad_points_1d.*0.5.*deltaXi,1,nquad_1d);
		quad_points(2,:) = quad_points(2,:) + repmat(eta1+0.5*deltaEta,1,nquad);

        % evaluate the geometry mapping and its derivatives
		
        F = nurb_eval(nurb,nurb.coords,2,quad_points);
		
        dF = nurb_derv_eval(nurb,nurb.coords,2,quad_points);
		
        jac_det = zeros(1,nquad);
		
        for i=1:nquad
            % evaluate the absolute value of the Jacobian's determinant
            jac_det(i) = abs(det(dF(:,:,i)));
        end

        % Compute the NURBS' area
		
        nurb_area = nurb_area + sum(quad_weight.*jac_det);


		% print nurb_area for testing purposes
		%disp("Nurb Area: " + nurb_area)

		%% Exercise 5
		
        % evaluate nurb basis functions and their derivatives
        coeffs = zeros(nshp, nurb.number(1), nurb.number(2));
        p = nurb.order(1) - 1;
        q = nurb.order(2) - 1;
        i0 = bspline_findspan(nurb.number(1), p, quad_points(1,:), nurb.knots{1});
        j0 = bspline_findspan(nurb.number(2), q, quad_points(2,:), nurb.knots{2});
        connectivity = zeros(ndof,nshp);
		
        % assemble coeffs and connectivity
        for i = 1:1:nurb.order(1)
            for j = 1:1:nurb.order(2)
                % Fill with ones, otherwise stays zero
                coeffs((j - 1) * (p + 1) + i, i0 - p + i - 1, j0 - q + j - 1) = 1;
                % For all degrees of freedom:
                for k=1:1:ndof
                    connectivity(k, (j - 1) * (p + 1) + i) ...
						= ndof * ((j0 - q + j - 2) * nurb.number(1) + i0 - p + i - 2) + k;
                end
            end
		end
		
		%disp(connectivity')

        if (use_bspline)
            % use B-Splines basis functions
            S = bspline_eval(nurb, coeffs, nshp, quad_points);
            % derivatives in the reference domain
            dS_ref = bspline_derv_eval(nurb, coeffs, nshp, quad_points);
        else
            % use NURBS basis functions

            % remember that nurb_(derv_)eval needs the coefficients
            % premultiplied with the weight
            % Remark: This step is only introduced to allow to test for
            % partition of unity of the NURBS basis functions. From an
            % algebraic point of view it does not matter whether we
            % premultiply the coefficients, it only introduces overhead in
            % the sense that the solution of the linear equation system
            % also needs to be mutiplied again with the weights before
            % plotting.
            for i=1:nshp
                coeffs(i,:,:) = coeffs(i,:,:) .* nurb.weights;
            end
            S = nurb_eval(nurb, coeffs, nshp, quad_points);
            % derivatives in the reference domain
            dS_ref =nurb_derv_eval(nurb, coeffs, nshp, quad_points);
        end
        dS_phys = zeros(nshp,2,nquad);

        % derivatives in the physical domain
        for i=1:nquad
            for j=1:nshp
                % compute dS_phys with the help of equation (3)
                % Use the '/' operator to compute the left inverse of the jacobian.
                dS_phys(j,:,i) = (dS_ref(j,:,i) / dF(:,:,i))';
            end
        end

        % check partition of unity
        testSum = sum(S,1);
        testSumDervRef = sum(dS_ref,1);
        testSumDervPhys = sum(dS_phys,1);
		for i=1:nquad
            assert(abs(testSum(1,i) - 1.0) < 10^-10);
            for j=1:2
                assert(abs(testSumDervRef(1,j,i)) < 10^-10);
				assert(abs(testSumDervPhys(1,j,i)) < 10^-10);
            end
		end
		
		%% Exercise 6
		for ishp = 1:nshp
			
			V = S(ishp,:);
			dV = squeeze(dS_phys(ishp,:,:));
			
			row = connectivity(:,ishp);
			
			for jshp = 1:nshp
				
				U = S(jshp,:);
				dU = squeeze(dS_phys(jshp,:,:));
				
				col = connectivity(:,jshp);
				
				% call element matrix
				matloc = element_matrix(ndof, quad_points, quad_weight, jac_det, F, U ,dU ,V ,dV);
				
				% add matloc to mat. Be aware of the connectivity array
				mat(row,col) = mat(row,col) + matloc;
			end
			
			% call rhs_function
			rhsloc = rhs_function(ndof, quad_points, quad_weight, jac_det, F, V, dV);
			% add rhsloc to rhs
			rhs(row) = rhs(row) + rhsloc;
		end
		
	end
	
end
end
