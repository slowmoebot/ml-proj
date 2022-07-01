clc
clear
rng(1)

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

addpath IGA\
load nets\net_min2x2_107rmse.mat
refinement_level = 5;

fig=figure(1);
clf(1)
n=3;
m=2;
tiledlayout(n,m,TileSpacing="tight",Padding="tight")

for i = 1:n

%% Create Random nurb

nurb = generate_random_nurb();
nurb = nurb_knot_refinement(nurb, refinement_level);

%% Prepare data

tic 
input(:,:,1) = squeeze(nurb.coords(1,:,:));
input(:,:,2) = squeeze(nurb.coords(2,:,:));
input(:,:,3) = squeeze(nurb.weights(1,:,:));

input(:,:,1) = ((input(:,:,1) ./ input(:,:,3))+0.5)./3;
input(:,:,2) = ((input(:,:,2) ./ input(:,:,3))+0.5)./3;
input(:,:,3) = input(:,:,3)/2;

res_cnn = zeros(1,3+refinement_level,3+refinement_level);
res_cnn(1,2:end-1,2:end-1) = predict(net,round(input*255))./255;

timings_cnn(i) =toc;

tic
res_iga=solve_diffusion_iga(nurb);
timings_iga(i) = toc;

nexttile
draw_nurb_surf(nurb, [15 15], @(xi,eta) nurb_eval(nurb, res_cnn, 1, [xi; eta]));
title("CNN")
xlim([-0.5 2.5])
ylim([-0.5 2.5])


nexttile
draw_nurb_surf(nurb, [15 15], @(xi,eta) nurb_eval(nurb, res_iga, 1, [xi; eta]));
title("IGA")
xlim([-0.5 2.5])
ylim([-0.5 2.5])


name="results_08x08";
width = 8;
height = 8;
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [width, height],"Renderer","painters");

%print(fig, sprintf("figs/%s.pdf", name), '-dpdf','-fillpage');

end