clearvars; close all;
rng(1);

N = 10000;
pd2 = makedist('Uniform');%,'sigma',5)
pd3 = makedist('Beta','a',0.5,'b',0.5);
pd4 = makedist('Beta','a',2,'b',2);
pd5 = makedist('Beta','a',1,'b',20);
pd6 = makedist('Normal','mu',0.5,'sigma',0.1); pd6 = truncate(pd6,0,1);
pd7 = makedist('Normal','mu',0.25,'sigma',0.05); pd7 = truncate(pd7,0,1);
pd8 = makedist('Normal','mu',0.75,'sigma',0.05); pd8 = truncate(pd8,0,1);
pd9 = makedist('Normal','mu',0.5,'sigma',0.0001); pd9 = truncate(pd9,0,1);

% Draw Samples  (Reordered in decreasing BSIEntropy)
d1 = linspace(0, 1, N); d1 = d1';   d2 = sort(pd2.random(N,1)); 
d4 = sort(pd3.random(N,1));         d3 = sort(pd4.random(N,1)); 
d8 = sort(pd5.random(N,1));         d7 = sort(pd6.random(N,1)); 
%Mixture of Gaussians
d6 = sort([pd7.random(3*N/4,1); pd8.random(N/4,1)]);
%% Spike and Slab
d5 = sort([pd2.random(3*N/4,1); 0.5*ones(N/4,1)]);
%d9 = sort([pd2.random(3*N/4,1);pd9.random(N/4,1)]);
%d10 = sort([pd2.random(N/4,1); 0.5*ones(3*N/4,1)]);

% MATLAB COLORS
colors = ["#0072BD","#D95319","#EDB120","#7E2F8E",'k',"#77AC30",...
    "#4DBEEE", "#A2142F",'r','g'];           
step = 250;
NN = N/step;
BSIE = zeros(NN,8);
for m=1:8
    d = eval(strcat('d',string(m)));
    mind = 0; maxd=1;    
    K = length(d);
    for k=step:step:K
        % histogram bins
        [p, edges] = histcounts(d,linspace(mind,maxd,k+1));
        p = p./sum(p);
        % quantile bins
        quants = linspace(0,1,k+1);
        V = [mind quantile(d,quants(1,2:end-1)) maxd];
        q = abs(diff(V))/(maxd-mind);          % abs for error in deltas  
        pq= (p+q)/2;
        % B-S measures
        ixq = ~q; ixp = ~p; ix = or(ixp,ixq);
        p_=p; p_(ix)=[]; q_=q; q_(ix)= [];        
        cnt = k/step;
        pq_ = pq;  pq_(ixq) = []; q(ixq)=[];
        BSIE(cnt,m) = dot(q,log2(q))-dot(q,log2(pq_));
        % remove empty bins for remaining measures
        ix = ~p;
        p(ix) = [];   pq(ix) = [];
        BSIE(cnt,m) = 1- (BSIE(cnt,m) + dot(p,log2(p))-dot(p,log2(pq)))/2;        
    end
    figure(1);
    plot(step:step:K, BSIE(:,m),'-o',	'Color', colors(m),'MarkerFaceColor',...
        colors(m),'MarkerSize',10); hold on
    set(gca,'fontsize',32) 

    figure(2);
    [f,x] = ecdf(d);
    if m==1
        plot(x,f,'--','Color', colors(m), 'Linewidth',4); hold on
        set(gca,'fontsize',32) 
    else
        plot(x,f,'-','Color', colors(m), 'Linewidth',2.5); hold on
        set(gca,'fontsize',32) 
    end
%     figure(3);
% 
%     if m~=5
%         scatter(0.1*m*ones(N,1), d, 400, 'o','filled','MarkerFaceColor', colors(m),...
%         'MarkerFaceAlpha',0.05); 
%     else
%         d = d(d~=0.5);
%         scatter(0.1*m*ones(length(d),1), d, 400, 'o','filled','MarkerFaceColor', colors(m),...
%             'MarkerFaceAlpha',0.005); 
%         scatter(0.1*m, 0.5, 400, 'o','filled','MarkerFaceColor', colors(m)); 
%     end
%     hold on
%     xlim([0 0.9])
%     xticks([])
%     ylim([0 1])
%     yticks([0.0 1.0])
%     ylabel('Data Sets')
%     set(gca,'fontsize',36) 
%     set(gca,'YAxisLocation','right')
 % Make the horizontal data distributions
    figure(4);
    if m~=5
        scatter(d, 0.1*m*ones(N,1), 400, 'o','filled','MarkerFaceColor', colors(m),...
        'MarkerFaceAlpha',0.05); 
    else
        d = d(d~=0.5);
        scatter(d, 0.1*m*ones(length(d),1), 400, 'o','filled','MarkerFaceColor', colors(m),...
            'MarkerFaceAlpha',0.005); 
        scatter(0.5, 0.1*m, 400, 'o','filled','MarkerFaceColor', colors(m)); 
    end
    hold on
    ylim([0 1])
    yticks([])
    xlim([0 1])
    xticks([0.0 1.0])
    xlabel('Data Sets')
    set(gca,'fontsize',32) 

end

figure(1); 
legend('$\mathcal{EP}$','$\mathcal{U}[0,1]$','$Beta(2,2)$','$Beta(\frac{1}{2},\frac{1}{2})$',...
    '$\mathcal{S}\mathcal{S}$','$\mathcal{G}\mathcal{M}$','$\mathcal{N}(\frac{1}{2},\frac{1}{100})$',...
    '$Beta(1,20)$','interpreter','latex','location','eastoutside','fontsize',32);
title({'Boltzmann-Shannon','Interaction Entropy'})
xlabel('Partition Size $K$','interpreter','latex');
ylabel('$1-JSD({\bf p}||{\bf q})$','interpreter','latex');
ylim([0 1])
yticks(0:0.5:1)
%axis('square')
hFig = figure(1);
set(hFig,'position', [100 100 1200 800]); 
% 
figure(2); 
legend('$\mathcal{EP}$','$\mathcal{U}[0,1]$','$Beta(2,2)$','$Beta(\frac{1}{2},\frac{1}{2})$',...
    '$\mathcal{S}\mathcal{S}$','$\mathcal{G}\mathcal{M}$','$\mathcal{N}(\frac{1}{2},\frac{1}{100})$',...
    '$Beta(1,20)$','interpreter','latex','location','eastoutside','fontsize',32);ylim([0 1])
title('Empirical CDFs')
ylabel('$Prob(X\leq x)$','interpreter','latex')
yticks([0 1])
xlabel('$x$','interpreter','latex')
hFig = figure(2);
set(hFig,'position', [100 100 1000 700]); 

% figure(3);
% 
% hFig = figure(3);
% set(hFig,'position', [1236,129,1000,225]); 

figure(4);

hFig = figure(4);
set(hFig,'position', [200,100,1500,450]); 
