function [unpruned,vstar,nfevals,coneFound,alpha,w,v,dimension,condition] = expandTreePruning(A,D,parents,nfevals,epsilon)

% A: Cell array with N matrices Ai
% D: Depth level in the binary tree (D = 1: Root)


N = length(A);
Aaux = cell2mat(A);
% Problem dimensions
r = length(A);
n = size(A{1},2);

nNodes = 2^D; % Number of nodes to be evaluated in each tree

coneFound = false; % flag to indicate wether a cone has been found
alpha = [];
w = [];
v = [];
condition = [];
dimension = [];
% W = [A1 A2 ... Ar]
% Li = [0 0 ... I ... 0]
% Matrix formed with r - 1 null blocks of dimension (n x n) and an (n x n)
% identity matrix at the ith block
auxZeros = zeros(n,n*r);
In = eye(n);
W = A{1};
L{1} = auxZeros; L{1}(1:n,1:n) = In;
for i = 2:r
    W = [W A{i}];
    L{i} = auxZeros; L{i}(1:n,1+(i-1)*n:i*n) = In;
end
F = null(W); % Matrix with columns forming a basis for the null space of W

unpruned = cell(n,1);
vstar = cell(n,1);

% For each variable (facet)
for i = 1:n
    % xi >= 1
    % Rays will be defined by points on the ith facet of the unit hybercube

    % Must choose a leaf in each of the n-1 binary trees at depth D
    % Nodes in each tree are numbered from 1 (root) to 2^(D+1)-1
    % At depth D, the nodes are 2^D, (2^D+1), ..., 2^(D+1)-1
    % Root: D = 0
    notpruned_ind = 0;
    children = zeros(2^(n-1)*numel(parents{i}),1);
    for j = 1:numel(parents{i})
        if ~isempty(parents{i})
            children((j-1)*2^(n-1)+1:j*2^(n-1)) = getChildren(parents{i}(j),n,D);
        else
            unpruned{i} = [];
        end
    end
    for ind_choice_nodes = 1:numel(children)
        choice_nodes = children(ind_choice_nodes);
        repBase = dec2basenum(choice_nodes,nNodes,n-1);
        for tree = 1:n-1
            node(tree) = repBase(tree) + 2^D; % 0 --> 2^D
            [xmin(tree),xmax(tree)] = intgen(node(tree));
        end

        % Generate guide points for the rays
        R = zeros(2^(n-1),n);
        for j = 1 : 2^(n-1)
            bminmax = dec2basenum(j,2,n-1); % 0 or 1
            R(j,i) = 1; % All points lie on the ith facet
            for tree = 1:n-1
                if bminmax(tree) == 0
                    xlim(tree) = xmin(tree);
                else
                    xlim(tree) = xmax(tree);
                end
            end
            R(j,[1:i-1 i+1:end]) = xlim;
        end

        [Haux,baux] = minCone(R);

        H = Haux;
        b = baux;

        brep = repmat(b,N,1);
        Hrep = H;
        for dummy = 2:N
            Hrep = blkdiag(Hrep,H);
        end

        % sum_{i=1}^{n} xi >= 1
        I = eye(n);
        ei = I(i,:);
        ei_times_sum = repmat(ei,1,N);
        Hrep = [Hrep; -ei_times_sum];
        brep = [brep; -1];
 
        nfevals = nfevals+1;
        [sStar,vmat] = testNullSpace(A,H,b,i,F,L);
        if sStar <= 0
            notpruned_ind = notpruned_ind + 1;
            unpruned{i} = [unpruned{i}, choice_nodes];
            vstar{i}{notpruned_ind} = vmat;
            v_vectors_as_columns = vmat; % matrix with the base vectors v_i as columns
            v = mean(v_vectors_as_columns,2); % mean of all columns
            alpha = v.'*v_vectors_as_columns; % projection of each v_i onto v to find alpha_i
            alpha = alpha/sum(alpha); % sum alpha_i = 1
            w = cell2mat(A)*reshape(kron(alpha,v),N*n,1); % w = sum_i(A_i * alpha_i * v)
            condition = norm(w,2)/norm(v,2); % calculate condition
            if condition < epsilon
               coneFound = true;
               dimension = i;
               return 
            end
        end
    end
end