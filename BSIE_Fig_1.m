clearvars; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%See circulargraph-master for graphs
%  Some mild changes were made to the original 
%   to highlight isolated nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = copper(40);
cmap = cmap(7:31,:);

k = 25;             % Number of bins used for histograms/quantiles

% Create Logistic Data
x = 0.1;
% for nn=1:24       % uncomment for subsample of N=25
for nn=1:9999
    x = [x; 4*x(end).*(1-x(end))];
end

quants = linspace(0,1,k+1);
V = [0 quantile(x,quants(1,2:end-1)) 1];        % Enforce end points [0,1]
q = diff(V);
[p,Eh] = histcounts(x,linspace(0,1,k+1));
p = p./sum(p);    




figure('units','centimeters','position',[30,15,20,10]);
bar(p, 'FaceColor', "#0072BD")
title({'Frequency-Based Distribution','(Histogram)'})
xlabel('Bins')
xticks([])
ylabel('{\bf p}')
%ylim([0, 0.55])
%yticks([0 0.1 0.2 0.3 0.4 0.5])
set(gca,'fontsize',20)

figure('units','centimeters','position',[30,15,20,10]);
bar(q,'FaceColor', "#EDB120")
title({'Measure-Based Distribution','(Quantiles)'})
xlabel('Quantiles')
xticks([])
ylabel('{\bf q}')
%ylim([0, 0.15])
%yticks([0 0.1])
set(gca,'fontsize',20)

figure;
ptr = scatter(x,0.01*ones(size(x)),75,'ok','MarkerFaceColor','k','MarkerEdgeColor','none');
ptr.MarkerFaceAlpha = 0.1;
hold on
for i=1:length(Eh)
    plot([Eh(i) Eh(i)], [-0.5 -0.1],'Color','#0072BD','LineWidth',1.5)
    plot([V(i) V(i)], [.1 0.5],'Color','#EDB120','LineWidth',1.5)
end
xlim([0 1])
xticks([0 1])
box off
%ylim([-0.5 0.5])
yticks([])
pbaspect([10 5 1])
legend({'$\rho(x)$','Data Sample','Histogram','Quantile'},'Orientation','horizontal','interpreter', 'latex');
%legend({'Data Sample','Quantile Bins'},'Orientation','horizontal','interpreter', 'latex');
set(gca,"FontSize",22)
p(~p)=[];
Hp = -dot(p,log(p))
q(~q)=[];
Hq = -dot(q,log(q))
