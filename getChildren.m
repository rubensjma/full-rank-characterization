function children = getChildren(parent,n,D)
%getChildren enumerates the children of a cone
% parent: index of the parent node
% n: dimension of the domain of the matrices
% D: depth in the three
% children: indices of the children

ndigits = n-1; % number of digits of the number that represents the position of the child in the (n-1) dimensional matrix
indicesParentString = dec2basenum(parent,2^(D-1),n-1); % number in base D-1
% shifts = dec2base([0:2^(n-1)-1].',2,ndigits); % all shifts of indices that result in other children
shifts = dec2bin([0:2^(n-1)-1].',ndigits); % all shifts of indices that result in other children
for i = 1:n-1
    digitsParent(1,i) = indicesParentString(n-i);
    digitsShifts(:,i) = str2num(shifts(:,n-i));
end

indicesChildren = repmat(2*digitsParent,2^(n-1),1) + digitsShifts; 

% Children integer equivalent from matrix indices
children = base2decnum(indicesChildren,2^D);

end