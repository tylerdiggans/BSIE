close all

rng(1);

pd2 = makedist('Uniform');%,'sigma',5)
pd3 = makedist('Beta','a',0.5,'b',0.5);
pd4 = makedist('Beta','a',2,'b',2);
pd5 = makedist('Beta','a',1,'b',20);
pd6 = makedist('Normal','mu',0.5,'sigma',0.1); pd6 = truncate(pd6,0,1);
pd7 = makedist('Normal','mu',0.25,'sigma',0.05); pd7 = truncate(pd7,0,1);
pd8 = makedist('Normal','mu',0.75,'sigma',0.05); pd8 = truncate(pd8,0,1);

% % Draw Samples
% N = 100000;
% D1 = linspace(0, 1, N);     D2 = sort(pd2.random(N,1)); 
% D3 = sort(pd3.random(N,1)); D4 = sort(pd4.random(N,1)); 
% D5 = sort(pd5.random(N,1)); D6 = sort(pd6.random(N,1));
% D7 = sort([pd7.random(3*N/4,1); pd8.random(N/4,1)]);
% 
% %% Spike and Slab
% D8 = sort([pd2.random(3*N/4,1);0.5*ones(N/4,1)]);
% % MATLAB COLORS
colors = ["#0072BD","#D95319","#EDB120","#7E2F8E",'k',"#77AC30",...
    "#4DBEEE", "#A2142F",'r','g'];                   

Inputs = [40, 100,500,1000,5000,10000,20000,50000,100000];

BSIE = zeros(length(Inputs),8);

for n=1:length(Inputs)
    % Draw Samples
    NN = Inputs(n);             
    d1 = linspace(0, 1, NN);     d2 = sort(pd2.random(NN,1)); 
    d4 = sort(pd3.random(NN,1)); d3 = sort(pd4.random(NN,1)); 
    d8 = sort(pd5.random(NN,1)); d7 = sort(pd6.random(NN,1));
    d6 = sort([pd7.random(3*NN/4,1); pd8.random(NN/4,1)]);
%% Spike and Slab
    d5 = sort([pd2.random(3*NN/4,1);0.5*ones(NN/4,1)]);
%     d1 = linspace(0,1, NN);      d2 = datasample(D2,NN);
%     d3 = datasample(D3,NN);       d4 = datasample(D4,NN);
%     d5 = datasample(D5,NN);       d6 = datasample(D6,NN);
%     d7 = datasample(D7,NN);       d8 = datasample(D8,NN);
    for m=1:8
        d = eval(strcat('d',string(m)));
        mind = 0; maxd=1;
        % histogram bins
        [p, edges] = histcounts(d,linspace(mind,maxd,NN+1));
        p = p./sum(p);
        % quantile bins
        quants = linspace(0,1,NN+1);
        V = [mind quantile(d,quants(1,2:end-1)) maxd];
        q = abs(diff(V))/(maxd-mind);          % abs for error in deltas
        pq = (p+q)/2;
        % B-S measures
        ixq = ~q;
        pq_ = pq;  pq_(ixq) = []; q(ixq)=[];
        BSIE(n,m) = dot(q,log2(q))-dot(q,log2(pq_));
        % remove empty bins for remaining measures
        ix = ~p;
        p(ix) = [];   pq(ix) = [];
        BSIE(n,m) = 1- (BSIE(n,m) + dot(p,log2(p))-dot(p,log2(pq)))/2;        
    end
end

for m=1:8
    figure(1);
    semilogx([40,100,500,1000,5000,10000,20000,50000,100000],BSIE(:,m),'-o','Color', colors(m),'MarkerFaceColor',...
        colors(m),'MarkerSize',10, 'Linewidth',2); hold on
end

figure(1);
title({'Boltzmann-Shannon','Interaction Entropy'})
xlabel('Sample size ($N$)','interpreter','latex')
ylabel('$BSIE(X_N)$','interpreter','latex')
legend('$\mathcal{EP}$','$\mathcal{U}[0,1]$','$Beta(2,2)$','$Beta(\frac{1}{2},\frac{1}{2})$',...
    '$\mathcal{S}\mathcal{S}$','$\mathcal{G}\mathcal{M}$','$\mathcal{N}(\frac{1}{2},\frac{1}{100})$',...
    '$Beta(1,20)$','interpreter','latex','location','eastoutside','fontsize',32);
ylim([0 1])
yticks([0 1])
xticks([100,1000,10000,100000])
set(gca,'fontsize',32) 
hFig = figure(1);
set(hFig,'position', [100 100 1200 800]); 