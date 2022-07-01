clear
clc

load nets\net_min2x2_107rmse.mat
layerg01 = layerGraph(net);
load nets\net_min1x1_107rmse.mat
layerg02 = layerGraph(net);
load nets\net23x23_min3x3_293rmse.mat
layerg03 = layerGraph(net);
load nets\net23x23_crop_296rmse.mat
layerg04 = layerGraph(net);


set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

fig=figure(1);
clf(1)
tiledlayout(1,4,TileSpacing="tight",Padding="tight")

nexttile
plot(layerg01)
title("Low refinement: A")
axis off

nexttile
plot(layerg02)
title("Low refinement: B")
axis off

nexttile
plot(layerg03)
title("High refinement: C")
axis off

nexttile
plot(layerg04)
title("High refinement: D")
axis off



name="networkgraphs";
width = 14.5;
height = 8;
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [width, height]);

print(fig, sprintf("figs/%s.pdf", name), '-dpdf','-fillpage');

