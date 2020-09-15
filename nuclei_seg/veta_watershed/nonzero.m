function y = nonzero(x)

switch class(x)
    case {'single'}
        epsilon = 1e-5;
    case {'double'}
        epsilon = 1e-200;
    otherwise
        error('X must be a floating point.');
end

y = x;
y(y==0) = epsilon;