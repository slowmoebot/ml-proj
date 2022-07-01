clc
clear
rng(1)

addpath IGA\

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig=figure(1);
clf(1)
n=3;
tiledlayout(n,3,"TileSpacing","tight","Padding","tight")

for i = 1:n
    nurb = generate_random_nurb();
    ref_nurb = nurb_knot_refinement(nurb,5);
    res_weights = solve_diffusion_iga(ref_nurb);


    nexttile

    draw_nurb_surf(ref_nurb, [30 30], @(xi,eta) nurb_eval(ref_nurb, res_weights, 1, [xi; eta]));
    colorbar off

    nexttile
    box on
    data.x = squeeze(ref_nurb.coords(1,:,:));
    data.y = squeeze(ref_nurb.coords(2,:,:)); 
    data.w = squeeze(ref_nurb.weights(1,:,:));

    data.sol = squeeze(res_weights(1,:,:));

    image(:,:,1)=((data.x ./ data.w)+0.5)./3;
    image(:,:,2)=((data.y ./ data.w)+0.5)./3;
    image(:,:,3)=data.w/2;

    imshow(image)

    nexttile
    
    imshow(data.sol(2:end-1,2:end-1));
    colormap("parula")
end

name="example_solutions";
width = 10;
height = 8;
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [width, height],"Renderer","painters");

print(fig, sprintf("figs/%s.pdf", name), '-dpdf','-fillpage');