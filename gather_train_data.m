clc
clear
addpath IGA\

n_data = 1000-515;
refinement_level = 5;

for i = 1:n_data

    nurb.coords = zeros(3,3,3);
    nurb.weights = zeros(1,3,3);
    nurb.coords(1,:,:) = reshape([0.00000  1.00000  2.00000  0.00000  1.00000  2.00000  0.00000  1.00000  2.00000], 3, 3);
    nurb.coords(2,:,:) = reshape([0.00000  0.00000  0.00000  1.00000  1.00000  1.00000  2.00000  2.00000  2.00000], 3, 3);
    nurb.coords(1:2,:,:) = nurb.coords(1:2,:,:) - 0.5 + rand(2,3,3);
    nurb.weights = rand(1,3,3) * 2;
    nurb.coords(1,:,:) = nurb.coords(1,:,:) .* nurb.weights;
    nurb.coords(2,:,:) = nurb.coords(2,:,:) .* nurb.weights;
    nurb.number = [ 3 3 ];
    nurb.order = [ 3 3 ];
    nurb.knots{1} = [0 0 0 1 1 1];
    nurb.knots{2} = [0 0 0 1 1 1];

    data(i).x_init = squeeze(nurb.coords(1,:,:));
    data(i).y_init = squeeze(nurb.coords(2,:,:)); 
    data(i).w_init = squeeze(nurb.weights(1,:,:));

    ref_nurb = nurb_knot_refinement(nurb,20);

    data(i).x = squeeze(ref_nurb.coords(1,:,:));
    data(i).y = squeeze(ref_nurb.coords(2,:,:)); 
    data(i).w = squeeze(ref_nurb.weights(1,:,:));

    res_weights = solve_diffusion_iga(ref_nurb);

    data(i).sol = squeeze(res_weights(1,:,:));

    clc
    fprintf("Iteration %d/%d complete.\n",i,n_data)
end