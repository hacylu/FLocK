function [h2] = elipLog(p2, sigma_x, sigma_y, theta)

% returns a rotational
%   Laplacian of Gaussian filter of size HSIZE with standard deviation
%   SIGMA (positive). HSIZE can be a vector specifying the number of rows
%   and columns in H or a scalar, in which case H is a square matrix.
%   The default HSIZE is [5 5], the default SIGMA is 0.5.

% p2=[15 15];
% sigma_x = 4;
% sigma_y = 4;
%theta = 
% Laplacian of Gaussian
% first calculate Gaussian
siz   = (p2-1)/2;

a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;

[x,y] = meshgrid(-siz(2):1:siz(2),-siz(1):1:siz(1));
arg = ( - (a*(x).^2 + 2*b*(x).*(y) + c*(y).^2)) ;

%      [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
%      arg   = -(x.*x + y.*y)/(2*std2);

h = exp(arg);
h(h<eps*max(h(:))) = 0;

sumh = sum(h(:));
if sumh ~= 0,
    h  = h/sumh;
end;
% now calculate Laplacian
%h1 = h.*(x.*x + y.*y - 2*std2)/(std2^2);
lapl = ((2*a*x+2*b*y).^2-2*a)+((2*c*y+2*b*x).^2-2*c);
h1 = h.*lapl;
h2 = h1 - sum(h1(:))/prod(p2); % make the filter sum to zero