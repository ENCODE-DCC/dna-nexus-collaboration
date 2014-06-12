function kval = generateKernel(type,bw,normflag,varargin)
% generates the kernel window function
% function kval = generateKernel(type,bw,normflag,varargin)
% type: kernel type
% bw: kernel bandwidth (window half width)
% normflag: if set to true kernel is normalized
% varargin{1}: tukeyRatio (taper/fixed)

assert((bw>0),'generateKernel:argError','kernel bandwidth must be > 0');
winLen = 2*ceil(bw)+1;
kval = linspace(-1,1,winLen);
switch(type)
    case {'unif','uniform','rect','rectangular'}        
        kval = ones(winLen,1);
    case {'tria','triangular'}
        kval = 1-abs(kval);
    case {'epan','epanechnikov'}
        kval = (3/4)*(1-kval.^2);
    case {'quart','quartic','biwe','biweight'}
        kval = (15/16)*(1-kval.^2).^2;
    case {'triw','triweight','tric','tricube'}
        kval = (35/32)*(1-kval.^2).^3;
    case {'cosi','cosine'}
        kval = (pi/4)*cos(pi*kval/2);
    case {'gaus','gaussian'}        
        kval = normpdf(linspace(-3,3,winLen),0,1);
    case {'tuke','tukey'}
        tukeyRatio = 0.25;
        if nargin > 3
            if ~isnan(varargin{1})
                tukeyRatio = varargin{1};
            end
        end
        kval = tukeywin(winLen,tukeyRatio);
    otherwise
        error('generateKernel:argError', ...
            'Supported kernel types: uniform,gaussian,epanechnikov,triweight,cosine,gaussian,tukey');
end
if normflag
    kval = kval/sum(kval); % normalize the kernel coefficients
else
    kval = kval/max(kval); % unnormalized kernel
end
kval = double(kval(:));
end