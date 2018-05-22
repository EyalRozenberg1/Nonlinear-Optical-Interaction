function u = thomas_algorithm(a, b, c, d)

% This function uses tridiagonal matrix algorigthm (Thomas algorithm)
% to solve the new vector u of heat equation
%
% Visa Suomi
% Turku University Hospital
% 
% November 2017

cc = NaN(size(d));
dd = NaN(size(d));
u = NaN(size(d));

cc(1) = c / b;
dd(1) = d(1) / b;

for ind = 2:length(d)
    cc(ind) = c / (b - a * cc(ind-1));
    dd(ind) = (d(ind) - a * dd(ind - 1)) / (b - a * cc(ind-1));
end

clear ind

u(end) = dd(end);

for ind = length(dd)-1:-1:1
    u(ind) = dd(ind) - cc(ind) * u(ind+1);
end

clear ind