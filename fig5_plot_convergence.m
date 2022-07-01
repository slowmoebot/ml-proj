clear
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load nets\net_min2x2_107rmse.mat
stat01 = net_stats.TrainingRMSE;
load nets\net_min1x1_107rmse.mat
stat02 = net_stats.TrainingRMSE;
load nets\net23x23_min3x3_293rmse.mat
stat03 = net_stats.TrainingRMSE;
load nets\net23x23_crop_296rmse.mat
stat04 = net_stats.TrainingRMSE;


fig=figure(1);
clf(1)
hold on
box on
grid on

plot(stat01,"-","LineWidth",2)
plot(stat02,"-","LineWidth",2)
plot(stat03,"-","LineWidth",2)
plot(stat04,"-","LineWidth",2)

legend(["Low refinement: A","Low refinement: B","High refinement: C","High refinement: D"])
xlabel("Iteration")
ylabel("RMSE Training Set")

name="convergence";
width = 10;
height = 6;
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [width, height]);

print(fig, sprintf("figs/%s.pdf", name), '-dpdf','-fillpage');