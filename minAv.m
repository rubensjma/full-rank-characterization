function [J,vstar] = minAv(Q,Hrep,brep,opts)

f = zeros(size(Q,1),1);
[vstar,J] = quadprog(Q,f,Hrep,brep,[],[],[],[],[],opts);
