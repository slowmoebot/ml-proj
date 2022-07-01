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
m=3;
tiledlayout(n,m,"TileSpacing","tight","Padding","tight")

for i = 1:n*m
    nexttile
    box on
    nurb = generate_random_nurb();
    draw_nurb_surf(nurb,[30 30], @(xi,eta) 1)
    colorbar off
end

name="example_domains";
width = 10;
height = 8;
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [width, height]);

print(fig, sprintf("figs/%s.pdf", name), '-dpdf','-fillpage');