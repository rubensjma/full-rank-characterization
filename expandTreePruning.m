function [Jmin,vstarmin,unpruned,nfevals] = expandTreePruning(A,D,parents,nfevals,tol_zero)

% A: Cell array with N matrices Ai
% D: Depth level in the binary tree (D = 1: Root)


N = length(A);
Aaux = cell2mat(A);
n = size(A{1},2);

Q = 2*Aaux'*Aaux;
Q = (Q + Q.')/2; % avoid numerical problems regarding symmetry

% Regularization for PD Q
evalues = eig(Q);
minevalue = min(evalues);
if minevalue < 0
    Q = Q - eye(size(Q,1))*minevalue;
end

nNodes = 2^D; % Number of nodes to be evaluated in each tree

unpruned = cell(n,1);

% For each variable (facet)
for i = 1:n
    % xi >= 1
    % Rays will be defined by points on the ith facet of the unit hybercube

    % Must choose a leaf in each of the n-1 binary trees at depth D
    % Nodes in each tree are numbered from 1 (root) to 2^(D+1)-1
    % At depth D, the nodes are 2^D, (2^D+1), ..., 2^(D+1)-1
    % Root: D = 0
    J = [];
    vstar = {};
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
        notpruned_ind = notpruned_ind + 1;
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
        opts = optimoptions('quadprog','Display','none');
        [J(notpruned_ind),vstar{notpruned_ind}] = minAv(Q,Hrep,brep,opts);
        if J(notpruned_ind) <= tol_zero % prune if cost > tolerance
            unpruned{i} = [unpruned{i}, choice_nodes];
        end

    end
    
% % Build map between children and parents
% [children,parents] = getChildParent(D,n-1);
% 
% nNodes = 2^D; % Number of nodes to be evaluated in each tree
% 
% % For each variable (facet)
% for i = 1:n
%     % xi >= 1
%     % Rays will be defined by points on the ith facet of the unit hybercube
% 
%     % Must choose a leaf in each of the n-1 binary trees at depth D
%     % Nodes in each tree are numbered from 1 (root) to 2^(D+1)-1
%     % At depth D, the nodes are 2^D, (2^D+1), ..., 2^(D+1)-1
%     % Root: D = 0
% 
%     children_vec = reshape(children,1,size(children,1)*size(children,2));
%     if numel(pruned{i,D-1}) > 0
%         unpruned = children(~ismember(parents,pruned{i,D-1})); % children of unpruned nodes
%     else
%         unpruned = children_vec;
%     end
%     children_vec(ismember(children_vec,unpruned)) = [];
%     pruned{i,D} = children_vec;
%     J = [];
%     vstar = {};
%     notpruned_ind = 0;
%     for ind_choice_nodes = 1:numel(unpruned)
%         choice_nodes = unpruned(ind_choice_nodes);
% 
%         notpruned_ind = notpruned_ind + 1;
%         repBase = dec2basenum(choice_nodes,nNodes,n-1);
%         for tree = 1:n-1
%             node(tree) = repBase(tree) + 2^D; % 0 --> 2^D
%             [xmin(tree),xmax(tree)] = intgen(node(tree));
%         end
% 
%         % Generate guide points for the rays
%         R = zeros(2^(n-1),n);
%         for j = 1 : 2^(n-1)
%             bminmax = dec2basenum(j,2,n-1); % 0 or 1
%             R(j,i) = 1; % All points lie on the ith facet
%             for tree = 1:n-1
%                 if bminmax(tree) == 0
%                     xlim(tree) = xmin(tree);
%                 else
%                     xlim(tree) = xmax(tree);
%                 end
%             end
%             R(j,[1:i-1 i+1:end]) = xlim;
%         end
% 
%         [Haux,baux] = minCone(R);
% 
%         H = Haux;
%         b = baux;
% 
%         brep = repmat(b,N,1);
%         Hrep = H;
%         for dummy = 2:N
%             Hrep = blkdiag(Hrep,H);
%         end
% 
%         % sum_{i=1}^{n} xi >= 1
%         I = eye(n);
%         ei = I(i,:);
%         ei_times_sum = repmat(ei,1,N);
%         Hrep = [Hrep; -ei_times_sum];
%         brep = [brep; -1];
% 
%         nfevals = nfevals+1;
%         opts = optimoptions('quadprog','Display','none');
%         [J(notpruned_ind),vstar{notpruned_ind}] = minAv(Q,Hrep,brep,opts);
%         if J(notpruned_ind) > tol_zero % prune if cost > tolerance
%             pruned{i,D} = [pruned{i,D}, choice_nodes];
%         end
% 
% 
%     end
    %     for choice_nodes = 0:( nNodes^(n-1) - 1)
    %         % find the parent of a node
    %         parent = parents(children==choice_nodes);
    %         if find(pruned{i,D-1} == parent)
    %             % parent pruned -> children pruned
    %             pruned{i,D} = [pruned{i,D},choice_nodes];
    %         else
    %             notpruned_ind = notpruned_ind + 1;
    %             repBase = dec2basenum(choice_nodes,nNodes,n-1);
    %             for tree = 1:n-1
    %                 node(tree) = repBase(tree) + 2^D; % 0 --> 2^D
    %                 [xmin(tree),xmax(tree)] = intgen(node(tree));
    %             end
    %
    %             % Generate guide points for the rays
    %             R = zeros(2^(n-1),n);
    %             for j = 1 : 2^(n-1)
    %                 bminmax = dec2basenum(j,2,n-1); % 0 or 1
    %                 R(j,i) = 1; % All points lie on the ith facet
    %                 for tree = 1:n-1
    %                     if bminmax(tree) == 0
    %                         xlim(tree) = xmin(tree);
    %                     else
    %                         xlim(tree) = xmax(tree);
    %                     end
    %                 end
    %                 R(j,[1:i-1 i+1:end]) = xlim;
    %             end
    %
    %             [Haux,baux] = minCone(R);
    %
    %             H = Haux;
    %             b = baux;
    %
    %             brep = repmat(b,N,1);
    %             Hrep = H;
    %             for dummy = 2:N
    %                 Hrep = blkdiag(Hrep,H);
    %             end
    %
    %             % sum_{i=1}^{n} xi >= 1
    %             I = eye(n);
    %             ei = I(i,:);
    %             ei_times_sum = repmat(ei,1,N);
    %             Hrep = [Hrep; -ei_times_sum];
    %             brep = [brep; -1];
    %
    %             nfevals = nfevals+1;
    %             opts = optimoptions('quadprog','Display','none');
    %             [J(notpruned_ind),vstar{notpruned_ind}] = minAv(Q,Hrep,brep,opts);
    %             if J(notpruned_ind) > tol_zero % prune if cost > tolerance
    %                 pruned{i,D} = [pruned{i,D}, choice_nodes];
    %             end
    %
    %         end
    %
    %     end

    if isempty(J) || sum(J <= tol_zero) == 0
        Jmin{i} = NaN;
        vstarmin{i} = NaN;
    else
        Jmin{i} = J(J <= tol_zero);
        vstarmin{i} = {vstar{J <= tol_zero}};
    end

end