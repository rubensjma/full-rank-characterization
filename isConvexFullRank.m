function [minCond,alphaMin,coneMin,nfevals,unpruned,vstar,v,dimension,w,Dfinish] = isConvexFullRank(A,maxD,epsilon)
%isConvexSingular Determines if there is a convex combination of given
%matrices that is singular.
%   A: cell of N vertices where each vertex is a n by n matrix
%   maxD: maximum depth of each tree in the splitting algorithm ( > 1)
%   tol_cond: numerical tolerance for the norm of the vector given
% by the mean v = mean(v_i)
% alpha_i = v_i' * v / v'*v
% w = sum_i(A_i * alpha_i * v)
% condition = norm(w)/norm(v)
%   minCond: value of the stopping condition
%   nfevals: number of LPs solved

N = length(A); % number of matrices
n = size(A{1},2); % dimension of the space
nfevals = 0; % number of LPs solved
D = 2; % initial depth for search ( D >= 2 )
unpruned = cell(n,maxD); % list of pruned nodes
for i = 1:n
       unprunedParents{i} = [0:2^((D-1)*(n-1))-1]; 
       unpruned{i,D} = unprunedParents{i};
end
condition = inf(n,1); % value of norm(w)/norm(v)
coneFound = false;
while D < maxD && prod(cellfun(@isempty,unprunedParents)) == 0 && coneFound == false
    % loop stops if either the maximal depth has been reached, all cones
    % are pruned in each direction or norm(w)/norm(v) <= tol_cond

    % solve the problem with prunning for depth D
    [unprunedParents,vstar,nfevals,coneFound,alpha,w,v,dimension,condition] = expandTreePruning(A,D,unprunedParents,nfevals,epsilon);
    unpruned(:,D) = unprunedParents; % save for plot
    D = D + 1; % increase depth
%     ind_evaluate = 1;
%     table_vstar = [];
%     condition = inf;
%     alpha = [];
%     w = [];
%     for j = 1:n % for all dimensions
%         if ~isempty(unprunedParents{j})
%             for k = 1:numel(unprunedParents{j}) % all cones where cost is considered zero
%                 v_vectors_as_columns = vstar{j}{k}; % matrix with the base vectors v_i as columns
%                 v = mean(v_vectors_as_columns,2); % mean of all columns
%                 alpha(ind_evaluate,:) = v.'*v_vectors_as_columns; % projection of each v_i onto v to find alpha_i
%                 alpha(ind_evaluate,:) = alpha(ind_evaluate,:)/sum(alpha(ind_evaluate,:)); % sum alpha_i = 1
%                 w(:,ind_evaluate) = cell2mat(A)*reshape(kron(alpha(ind_evaluate,:),v),N*n,1); % w = sum_i(A_i * alpha_i * v)
%                 condition(ind_evaluate) = norm(w(:,ind_evaluate),2)/norm(v,2); % calculate condition
%                 ind_evaluate = ind_evaluate + 1;
%                 table_vstar = [table_vstar; j k];
%             end
%         end
%     end
end

Dfinish = D-1;
if condition < epsilon
    minCond = condition;
    alphaMin = alpha;
    vstar = vstar{dimension}{end};
    auxConeMin = unpruned{dimension,Dfinish};
    coneMin = auxConeMin(end);
else
    minCond = [];
    alphaMin = [];
    vstar = [];
    coneMin = [];
end
% if sum(condition < epsilon) >= 1
%     [minCond,minCondInd] = min(condition);
%     alphaMin = alpha(minCondInd,:);
%     vstar = vstar{table_vstar(minCondInd,1)}{table_vstar(minCondInd,2)}; % square matrix with the base vectors v_i as columns
%     auxConeMin = unpruned{table_vstar(minCondInd,1),Dfinish};
%     coneMin = auxConeMin(table_vstar(minCondInd,2));
% else
%     minCond = [];
%     alphaMin = [];
%     vstar = [];
%     coneMin = [];
% end
end