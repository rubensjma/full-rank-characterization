% This work was developed by 
% Marcelo Carvalho Minhoto Teixeira, marcelo.minhoto@unesp.br
% Roberto Kawakami Harrop Galvão, kawakami@ita.br
% Edvaldo Assunção, edvaldo.assuncao@unesp.br
% Rubens Junqueira Magalhães Afonso, rubensjm@ita.br

clear, clc, close all

% Uncomment for the desired example and run

%% Example 1, there is no singular combination

% example = 1;
% A{1} = [4     4    -6;
%       0    -2     2;
%      -4    -2    -5];
% A{2} = [1    -3    -3;
%      -10   -4    -1;
%      -4     6    -1];
% A{3} = [1     8     1;
%       0    -3    -2;
%      -5    -6    -2];

%% Example 2, inconclusive by the state-of-the-art. No combination exists

% example = 2;
% A{1} = [2	-2	2	0
% -1	-6	0	-1
% -4	-1	-1	-4
% -5	0	-6	5];
% A{2} = [-3	7	-1	-1
% -3	1	2	0
% -1	-3	0	3
% 0	1	1	2];
% A{3} = [0	2	7	3
% 0	0	1	-1
% 0	0	-4	-1
% -3	-1	2	0];


%% Example 3, non-square matrices, not full-row rank (consider alpha_i = 1/3, i = 1,2,3, as an example)

example = 3;
A{1} = [2     1    -2
     1     0     0
     2     1     2
     0     0     0];

A{2} = [-1     0     1
    -1    -1    -2
     0    -2     1
     2     1    -1];

A{3} = [-1    -1     3
     0     1     2
    -2     1    -3
    -3    -2    -2];

%% Parameters
maxD = 20; % stop if depth exceeds this value
epsilon = 1e-5; % numerical tolerance for the norm of the vector given
% by the mean v = mean(v_i)
% alpha_i = v_i' * v / v'*v
% w = sum_i(A_i * alpha_i * v)
% condition = norm(w)/norm(v)

n = size(A{1},2); % Dimension

%% Main routine
% If all components of Jmin = NaN, then all dimensions have been pruned and
% there is no convex combination
[minCond,alphaMin,coneMin,nfevals,unpruned,vstar,v,dimension,w,Dfinish] = isConvexFullRank(A,maxD,epsilon)
save(['example',int2str(example)]);