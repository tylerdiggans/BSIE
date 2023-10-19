clearvars; close all;
rng(1);

N1 = 100000; 
pd1 = makedist('Normal'); 
pd2 = truncate(pd1,-5,5);
pd1 = truncate(pd1,-3,3);
N2 = 1000; 
N3 = 200; 
pd3 = makedist('Uniform','lower',0,'upper',1);
N4 = 2000; 
pd4 = makedist('Generalized Pareto');  pd4 = truncate(pd4,1,1000);
N5 = 100; 

% Draw Samples
d1 = sort(pd1.random(N1,1)); 
d2 = sort(2*pd1.random(N2,1)+1);  % linear transform and different N
d3 = sort(pd2.random(N2,1)); 
d4 = sort(pd3.random(N3,1)); 
d5 = sort(pd4.random(N4,1)); 
d6 = sort(pd4.random(N5,1));
Ns = [N1,N2,N2,N3,N4, N5];
alphas = [0.01 0.2 0.2 0.2 0.1 0.1];
% MATLAB COLORS
colors = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30",...
    "#4DBEEE", "#A2142F",'k'];           
mind = [-3, -5, -5,0, 1, 1]; 
maxd = [3, 7, 5,1, 1000, 1000]; 
BSRE = zeros(6,1); BSGA = zeros(6,1); BSIE = zeros(6,1);
DSE = zeros(6,1); GPE = zeros(6,1);

for m=1:6
    d = eval(strcat('d',string(m)));
    % histogram bins
    [p, edges] = histcounts(d,linspace(mind(m),maxd(m),Ns(m)+2));
    p1 = p./sum(p);
    q1 = diff(edges)/(maxd(m)-mind(m));
    q1 = q1/sum(q1);

    % approx quantile bins
    V = [mind(m) d' maxd(m)];
    q2 = abs(diff(V))/(maxd(m)-mind(m));          % abs for error in deltas
    p2 = histcounts(d,V);
    p2 = p2./sum(p2);    
    pq = (p1+q2)/2;
    % B-S measures
    ixq = ~q2; ixp = ~p2; ix = or(ixp,ixq);
    p2_=p2; p2_(ix)=[]; q2_=q2; q2_(ix)= [];        
    BSGA(m) = dot(q2_,log(q2_))-dot(q2_,log(p2_));
    pq_ = pq;  pq_(ixq) = []; q2(ixq)=[];
    BSIE(m) = dot(q2,log(q2))-dot(q2,log(pq_));
    % remove empty bins for remaining measures
    ix = ~p1;
    p1(ix) = [];    q1(ix) = [];    pq(ix) = [];
    BSIE(m) = 1- (BSIE(m) + dot(p1,log(p1))-dot(p1,log(pq)))/(2*log(2));        
    BSRE(m) = dot(p1,log(p1))-dot(p1,log(q1));

    [p,edges] = histcounts(d,'BinMethod','scott');
    nn = length(edges)
    p = p/sum(p);
    p(~p)=[];
    DSE(m) = -dot(p,log(p));

    % quantile bins
    quants = linspace(0,1,nn);
    V = [mind(m) quantile(d,quants(1,2:end-1)) maxd(m)];
    q = abs(diff(V))/(maxd(m)-mind(m));          % abs for error in deltas
    q(~q)=[];
    GPE(m) = -dot(q,log(q));

    figure(1);
    [f,x] = ecdf(d);
    if m==1 || m==5
        plot((x-mind(m))/(maxd(m)-mind(m)),f,'--','Color', colors(m),...
            'Linewidth',4); hold on
    else
        plot((x-mind(m))/(maxd(m)-mind(m)),f,'-','Color', colors(m),...
            'Linewidth',2); hold on
    end
    set(gca,'fontsize',20) 

 % Make the horizontal data distributions
    figure(2);
    scatter(0.1*m*ones(Ns(m),1), (d-mind(m))/(maxd(m)-mind(m)), 400, 'o','filled','MarkerFaceColor', ...
                colors(m),'MarkerFaceAlpha',alphas(m)); 
    hold on
    xticks([])
    ylim([0 1])
    yticks([0.0 1.0])
    ylabel('Data Sets')
    set(gca,'fontsize',36) 
    set(gca,'YAxisLocation','right')
end
format long
BSGA
BSRE
BSIE
DSE
GPE
figure(1); 
hold on
plot(linspace(0,1,50),[0 linspace(0,1,49)],'--k')
%legend('$100000-\tau\mathcal{N}(0,1;-3,3)$',...
%        '$1000-2\left[\tau\mathcal{N}(0,1;-3,3)\right]+1$',...
%        '$1000-\tau\mathcal{N}(0,1;-5,5)$','$200-\mathcal{U}[0,1]$',...
%        '$2000-\tau\mathcal{PL}(2;1,1000)$',...
%        '$100-\tau\mathcal{PL}(2;1,1000)$',...
%        'location','eastoutside','interpreter','latex','fontsize',24);
ylim([0 1])
title('Empirical CDFs','interpreter','latex')
ylabel('$Prob(X\leq x)$','interpreter','latex')
yticks([0 1])
xlabel('$x$','interpreter','latex')
hFig = figure(1);
set(hFig,'position', [500 100 700 700]); 

figure(2);
xlim([0 0.7])
xticks([])
xticks([0.1 0.2 0.3 0.4 0.5 0.6])
xticklabels({'1','2','3','4','5','6'})
yticks([0 1])
yticklabels({'a','b'})
hFig = figure(2);
set(hFig,'position', [1236,129,225,1000]); 




