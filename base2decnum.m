function decimal = base2decnum(num,base)
%base2decnum Converts from base to decimal
%   num is a matrix, each column is a digit, each row is a number
ndigits = size(num,2);
factor = base.^[0:ndigits-1].';
decimal = num*factor;
end

