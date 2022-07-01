function [sol_coeffs] = solve_diffusion(nurb)

addpath IGA\

% degrees of freedom per shape function
ndof = 1;
[mat, rhs] = assemble_matrix_2d(false, ndof, nurb, @blk_dudv, @rhs_testfun);

n = nurb.number(1);
m = nurb.number(2);

all_dofs = 1:n*m;

% equation numbers that correspond to weighting functions with support on
% the boundary with eta=0
boundary = 1:n;
% boundary with eta=1
boundary = union(boundary,(m-1)*n:m*n);
% boundary with xi=0
boundary = union(boundary,1:n:m*n);
% boundary with xi=1
boundary = union(boundary,n:n:m*n);

% subtract the boundary set from all degrees
inner_dofs = setdiff(all_dofs,boundary);

% initialize solution vector to zero, including the boundary terms
solution = zeros(n*m,1);

% solve the linear system only for the inner degrees
solution(inner_dofs) = mat(inner_dofs,inner_dofs) \ rhs(inner_dofs);

% reshape the solution into a coefficients array
sol_coeffs = reshape(solution,ndof,n,m);
% remember that the nurb_eval routine needs the coefficients premultiplied
% with the weights
for i=1:ndof
    sol_coeffs(i,:,:) = sol_coeffs(i,:,:).*nurb.weights;
end


end