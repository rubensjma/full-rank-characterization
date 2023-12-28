function [sStar,vmat] = testNullSpace(A,H,b,m,F,L)

% Inputs:
% A : Cell array containing matrices A1, A2, ..., Ar
% H,b: H-rep representation of the cone C(J_m,n(d)) = {v : H v <= b}
% m : Facet of the unit hypercube under consideration
% F: Matrix with columns forming a basis for the null space of W = [A1 A2 ... Ar]
% L{i} = [0 0 ... I ... 0]: Matrix formed with r - 1 null blocks of dimension (n x n) and an (n x n) identity matrix at the ith block
%
% Output:
% The problem under consideration is feasible iff sStar <= 0
% vcell: Cell array comprising
% vcell{1} = v1Star, vcell{2} = v2Star, ..., vcell{r} = vrStar

% Problem dimensions
r = length(A);
n = size(A{1},2);
 
% % W = [A1 A2 ... Ar]
% % Li = [0 0 ... I ... 0]
% % Matrix formed with r - 1 null blocks of dimension (n x n) and an (n x n)
% % identity matrix at the ith block
% auxZeros = zeros(n,n*r);
% In = eye(n);
% W = A{1};
% L{1} = auxZeros; L{1}(1:n,1:n) = In;
% for i = 2:r
%     W = [W A{i}];
%     L{i} = auxZeros; L{i}(1:n,1+(i-1)*n:i*n) = In;
% end
% F = null(W); % Matrix with columns forming a basis for the null space of W

% The problem consists in determining whether there exists z s.t.
% H*(Li*F*z) <= b, i = 1, 2, ..., r
% -em'*(L1 + L2 + ... Lr)*F*z <= -1

nz = size(F,2); % Dimension of the z vector
Hz = H*L{1}*F;
bz = b;
for i = 2:r
    Hz = [Hz;H*L{i}*F];
    bz = [bz;b];
end

em = zeros(n,1); em(m) = 1;
sumL = L{1};
for i = 2:r
    sumL = sumL + L{i};
end
Hz = [Hz;-em'*sumL*F];
bz = [bz;-1];
% Number of constraints on z
mz = size(Hz,1);

% linprog formulation
% Vector of decision variables: [z;s]
% Cost: s
% Constraints: 
% Hz*z - 1_mz*s <= bz
%            -s <= 1

flp = zeros(nz+1,1); flp(end) = 1;
Alp = [Hz -ones(mz,1);
       zeros(1,nz) -1];
blp = [bz;1];

% % No initial solution
% opts = optimoptions('linprog','display','none');
% xlp = linprog(flp,Alp,blp,[],[],[],[],[],opts);

% Feasible initial solution with Gurobi
x0 = [zeros(mz,1); max([-1 -bz.'])]; % initializing with feasible solution
% model.A = sparse(Alp);
% model.B = blp;
% model.f = flp;
% model.x0 = x0;
% 
% xlp = gurobi(model)

lb = [];
ub = [];
ctype = repmat('U',mz+1,1); % ineqs <=
vartype = repmat('C',nz+1,1); % continuous opt variables 
sense = 1; % minimization
param.msglev = 1; % error messages only
[xlp, fmin, status, extra] = glpk (flp, Alp, blp, lb, ub, ctype, vartype, sense, param);

zStar = xlp(1:end-1);
sStar = xlp(end);
for i = 1:r
    vcell{i} = L{i}*F*zStar;
end
vmat = cell2mat(vcell);