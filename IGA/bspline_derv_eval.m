function dS = bspline_derv_eval(bspline, coeffs, dimS, points)
% This routine calculates for every points(:,i) the corresponding Jacobian of the
% B-Spline function and stores it in dS(:,:,i).

% bspline:                 holds the geometrical B-Spline structure
% coeffs(dimS,n1[,n2]):    the coefficients we should use for evaluation, which are not
%                          necessarily the control points
% dimS:                    dimension of the space the curve/surface/DOF is embedded in,
%                          which implies that the first index of coeffs should have at
%                          least that many entries.
% points(dim,size_points): the evaluation points with points(1,i) and points(2,i)
%                          being the $\xi$- and $\eta$-direction of the
%                          i-th evaluation point

dim = numel(bspline.order);
% Sanity checks
assert(dim == numel(bspline.number));
assert(size(points,1) == dim);
for i=1:dim
	assert(size(coeffs,i+1) == bspline.number(i))
end
assert(size(coeffs,1) >= dimS);

size_points = size(points,2);

% Preallocate dS
dS = zeros(dimS, dim, size_points);

if (dim == 1)
	% new data strucutre dbspline
	p = bspline.order(1)-1;
	dbspline.knots{1} = bspline.knots{1}(2:end-1);
	dbspline.order = bspline.order(1)-1;
	
	dbspline.number = bspline.number(1)-1;
	%dbspline.coords = (bspline.coords(:,2:end)-bspline.coords(:,1:end-1)).*p./(bspline.knots{1}(p+2:end-1)-bspline.knots{1}(2:end-p-1));
	dcoeffs = (coeffs(:,2:end)-coeffs(:,1:end-1)).*p./(bspline.knots{1}(p+2:end-1)-bspline.knots{1}(2:end-p-1));
	% evaluate
	dS(:,1,:)=bspline_eval(dbspline,dcoeffs,dimS,points);
	
	
elseif (dim == 2)
    % Compute derivative in $\xi$-direction:
    % Fill the dbspline structure as in the dim == 1 case, leaving the $\eta$-related
    % information unchanged. Store the result of bspline_eval in the first column
    % of dS.
    dbspline.number = [ bspline.number(1)-1 bspline.number(2)];
    dbspline.order = [ bspline.order(1)-1 bspline.order(2)];
    dbspline.knots{1}=bspline.knots{1}(2:end-1);
    dbspline.knots{2}=bspline.knots{2};
    p=bspline.order(1)-1;
    %dbspline.coords=zeros(dimS,bspline.number(1)-1,bspline.number(2));
    %dbspline.coords(:,:,:) = (bspline.coords(1:dimS,2:end,:)-bspline.coords(1:dimS,1:end-1,:)).*p./(bspline.knots{1}(p+2:end-1)-bspline.knots{1}(2:end-p-1));    

    
    dcoeffs1=zeros(dimS,bspline.number(1)-1,bspline.number(2));
    dcoeffs1(:,:,:) = (coeffs(1:dimS,2:end,:)-coeffs(1:dimS,1:end-1,:)).*p./(bspline.knots{1}(p+2:end-1)-bspline.knots{1}(2:end-p-1));
    dS(:,1,:)=bspline_eval(dbspline,dcoeffs1,dimS,points);

    % Compute derivative in $\eta$-direction:

    dbspline.number = [ bspline.number(1) bspline.number(2)-1];
    dbspline.order = [ bspline.order(1) bspline.order(2)-1];
    dbspline.knots{1}=bspline.knots{1};
    dbspline.knots{2}=bspline.knots{2}(2:end-1);
    q=bspline.order(2)-1;
    %dbspline.coords=zeros(dimS,bspline.number(1),bspline.number(2)-1);
    %dbspline.coords(:,:,:) = (bspline.coords(1:dimS,:,2:end)-bspline.coords(1:dimS,:,1:end-1)).*q./(bspline.knots{2}(q+2:end-1)-bspline.knots{2}(2:end-q-1));
    dcoeffs2=zeros(dimS,bspline.number(1),bspline.number(2)-1);
    denom = reshape(bspline.knots{2}(q+2:end-1)-bspline.knots{2}(2:end-q-1),1,1,[]);
    dcoeffs2(:,:,:) = (coeffs(1:dimS,:,2:end)-coeffs(1:dimS,:,1:end-1)).*q./denom;
    dS(:,2,:)=bspline_eval(dbspline,dcoeffs2,dimS,points);
end
end
