clc
clear

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

n_points = 101;
t= linspace(0,6*pi,n_points);
f= sin(t)+0.1*randn(1,n_points);

D = [t;f]';

fig = figure(1);
clf(1)
hold on
grid on
box on
n_cps = [4,6,8,10];
for n_cp = n_cps
    nurb = spline_fit(D,2,n_cp);
    
    points = linspace(0,1,101);
    
    evals = nurb_eval(nurb,nurb.coords,2,points);
    
    plot(evals(1,:),evals(2,:),"DisplayName",sprintf("$n=%d$",n_cp))
end
scatter(D(:,1),D(:,2),"DisplayName","Points")
legend("Location","northeast")
xlabel("$x$")
ylabel("$y$")


name="spline_regression";
width = 10;
height = 8;
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [width, height]);

print(fig, sprintf("figs/%s.pdf", name), '-dpdf','-fillpage');