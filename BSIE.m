function D = BSIE(x,limits, k)
x = sort(x);
if nargin<2
    %%% Use the extreme values as bounds (with default k)
    q = (x(2:end)-x(1:end-1))/(max(x)-min(x));
    k = length(x)-1;
elseif nargin<3
    %%% Use default k value, but add limits that are provided before
    %   computing quantiles
    if ~isnan(limits(1))
        if x(1)<limits(1)
            fprintf('ERROR: data outside limits')
            return
        else
            y = [limits(1); x];
        end
    else
        y = x;      %  if bound is nan skip
    end
    if ~isnan(limits(2)) 
        if x(end)>limits(2)
            fprintf('ERROR: data outside limits')
            return
        else
            y = [y; limits(2)];
        end
    elseif isnan(limits(1))
        y = x; %  if bound is nan skip
    end
    q = (y(2:end)-y(1:end-1))/(max(y)-min(y));
    k = length(y)-1;    
else
    %%% Geometry optional with k bins
    t = linspace(0,1,k+1)';
    s = quantile(x,t);
    l = s(2:end)-s(1:end-1);
    q = l./sum(l);    
end

%%%Histogram
t = linspace(min(x),max(x),k+1)';
s = discretize(x,t);
h = accumarray(s,1,[k 1]);
p = h./sum(h);

%%% Divergence
m = (p+q)./2;
ix = ~p; p(ix) = []; m_ = m; m_(ix) = [];
D = dot(p,log(p./m_));
ix = ~q; q(ix) = []; m(ix) = [];
D = 1-0.5*(D+dot(q,log(q./m)))/log(2);
end