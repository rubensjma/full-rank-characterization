function [xmin,xmax] = intgen(node)

% node = 1, 2, ...
% 1 = Root

xmin = -1;
xmax = 1;

aux = dec2bin(node);
t = aux(2:end);
L = length(t);
for k = 1:L
    xmid = (xmin + xmax)/2;
    if t(k) == '0'
        xmax = xmid;
    else
        xmin = xmid;
    end
end