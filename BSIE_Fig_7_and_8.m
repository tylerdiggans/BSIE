clearvars; close all; clc;

% This create the plots for Figures 6 and 7 in
% 
% There are 6 regression examples: 3 curves sampled using 2 distributions 
%  for x N samples are drawn from the distributions for x and the 
%  corresponding y values are computed.  Gaussian noise $/nu$ is added to 
%  the outputs (y's) to simulate measurement error

N = 10000;                  % Total Sample Size
NoiseSTD = 0.025;           % Amplitude of Gaussian Noise
n = 100;                    % Subsample Size

% Distributions for independent variable (\subseteq[-1,1] for consistency)
xExp = {'1/5*randn(N,1)','2*betarnd(4,0.5,[N,1])-1'};
xLatex = {'$X\subseteq\mathcal{N}\left(0,\frac{1}{25}\right)$',...
            '$X\subseteq 2 \cdot beta\left(4,\frac{1}{2}\right)-1$'};
xName = {'Normal','Beta'};

% Example function for Regression
funcExp = {'x.^3+x.^2','x.^8+4*x.^4-x.^2+x-2','sin(2*pi*x)'};
funcLatex = {'$y=x^3+x^2+\nu$','$y=x^8+4x^4-x^2+x-2+\nu $',...
                '$y=\sin\left(2\pi x\right)+\nu$'};
funcName = {'Cubic','Octic','Sine'};
Beta_True = [0, 0, 1, 1, 0, 0 , 0, 0, 0; ...
         -2, 1, -1, 0, 4, 0 , 0, 0, 1; ...
         0, 2*pi, 0, -41.3417, 0, 81.605249, 0, -76.70586, 0]';
            % Note: the interval of convergence for Sine is (-1,1)

% Set up linspace for plotting smooth polynomial fit extending past min/max 
xx = linspace(-1,1,100)';    
C = [ones(100,1), xx, xx.^2, xx.^3, xx.^4, ...
            xx.^5, xx.^6, xx.^7,xx.^8];

% Choose a random sample of the independent variable
for i = 1:2
    x = sort(eval(string(xExp(i)))); 
    xmin = x(1);    xmax = x(end);
    A = [ones(size(x)), x, x.^2, x.^3, x.^4, ...
                        x.^5, x.^6, x.^7,x.^8];
    BSI_0 = BSIE(x)
    for j=1:3
        func = eval(strcat('@(x)',string(funcExp(j))));
        % Measure dependent variable with noise
        y = func(x) + NoiseSTD*randn(N,1);  % \nu is gaussian noise ~0.1
        Beta_Total = pinv(A)*y;             % least squares fit all data

        % Plot Data, true underlying curve, and polynomial fit for data
        figure;
        scatter(x,y,'or');
        hold on;
        plot(xx,func(xx),'-b','LineWidth',1.5);
        yy = C*Beta_Total;        % Plot polynomial fit
        plot(xx,yy,'--m','LineWidth',2)
        xlim([-1,1]);
        %title({string(funcLatex(j)),...
        %    string(xLatex(i))},'fontsize',20,'Interpreter','latex')
        xlabel('x','Interpreter','latex')
        ylabel('y','Interpreter','latex')
        legend('Data','True','Poly($\beta$)','Interpreter','latex')
        set(gca,"FontSize",20)
        set(gcf,'Position',[100 100 750 500])
        saveas(gcf,strcat('../subsample/',string(funcName(j)),string(xName(i)),'.png'))

        for r=1:50000
            ix = sort(randperm(N,n));
            xsub = x(ix);
            ysub = y(ix);
            BSI_r(r,1) = BSIE(xsub);
            AA = [ones(n,1), xsub, xsub.^2, xsub.^3, xsub.^4, ...
                            xsub.^5, xsub.^6, xsub.^7,xsub.^8];
            Beta = pinv(AA)*ysub;
            E(r,1) = norm(Beta_True-Beta,1)/norm(Beta_True,1);
        end
        RE  = abs(BSI_r-BSI_0)/abs(BSI_0);
        t = linspace(0,max(RE),11)';
        Bins = discretize(RE,t);
        for s=1:max(Bins)
            ix = Bins==s;
            MRE(s,1) = mean(E(ix));
        end

        figure;
        hold on
        subplot(1,2,1)
        scatter(RE,E,10,'blue','filled', 'MarkerFaceAlpha', 0.25)
        hold on
        grid on
        tt = 0.5*(t(2:end)+t(1:end-1));         % midpoints for markers
        scatter(tt,MRE,'dr','LineWidth',1.5)
        % fit an exponential curve to the stable (non-outlier) values
        f = fit(tt(1:end-2),MRE(1:end-2),'exp1');
        plot(f,tt,MRE,'dr')
        hFig=findall(0,'type','figure');
        hLeg=findobj(hFig(1,1),'type','legend');
        set(hLeg,'visible','off')
        xlabel('$RE\left(BSIE(X^i)\right)$','Interpreter','latex')
        ylabel('$RE\left(\hat{\beta}^i\right)$','Interpreter','latex')
        set(gca,"FontSize",20)

        subplot(1,2,2)
        scatter(tt,MRE,'dr','LineWidth',1.5);
        hold on
        grid on
        plot(f,tt,MRE,'dr')
        hFig=findall(0,'type','figure');
        hLeg=findobj(hFig(1,1),'type','legend');
        set(hLeg,'visible','off')
        xlabel('$RE\left(BSIE(X^i)\right)$','Interpreter','latex')
        ylabel('$RE\left(\hat{\beta}^i\right)$','Interpreter','latex')
        set(gca,"FontSize",20)
        %sgtitle({strcat(string(funcName(j)),' Data Set'),...
        %    string(xLatex(i))},'fontsize',20,'Interpreter','latex')
        set(gcf,'Position',[100 100 750 500])
        saveas(gcf,strcat('../subsample/',string(funcName(j)),string(xName(i)),'_BSI.png'))
    end
end