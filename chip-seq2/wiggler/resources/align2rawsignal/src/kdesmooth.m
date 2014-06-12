function [x] = kdesmooth(x,halfWidth,kernelType,normFlag,varargin)
% Smoothes x using kernel density smoothing
% function [x] = kdesmooth(x,bandWidth,kernelType,normFlag,varargin)
% x: input vector
% halfWidth: half the kernel window size
% kernelType: name of smoothing kernel
% normFlag: if set to true kernel smoothing is equivalent to smooth averaging. If not set, it is
%           equivalent to a smooth sum.
% varargin{1}: tukeyRatio

sizeX = size(x);
assert((min(sizeX)==1),'input vector x cannot be a matrix');

isRowVector = (sizeX(1)==1);

% Convert x to column vector
if isRowVector
    x = x(:);
end

% Check for tukeyRatio
tukeyRatio = NaN;
if nargin > 4
    tukeyRatio = varargin{1};
end

% Generate kernel
kval = generateKernel(kernelType,halfWidth,normFlag,tukeyRatio);
normKval = generateKernel(kernelType,halfWidth,1,tukeyRatio);

% Compute convolution
denom = isfinite(x); % Check for NaNs or Infs in x
hasNaN = any(~denom);
if hasNaN
    x(~denom) = 0; % Replace NaN and Infs with 0
end
% The division normalizes for missing values represented by NaNs or Infs
x = conv2( x , kval , 'same' ) ./ conv2( double(denom) , normKval , 'same' );

% Reset original NaN values to NaN
x(~denom) = NaN;

if isRowVector
    x = x';
end

end