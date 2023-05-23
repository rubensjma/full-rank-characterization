clear, clc, close all
disp('Requires function remark3.m available in Appendix 2 of the following paper:')
disp('Comments on `Less conservative conditions for robust LQR-state-derivative controller design: an LMI approach’ and new sufficient LMI conditions for invertibility of a convex combination of matrices.');
disp('https://doi.org/10.1080/00207721.2021.2023689')

A1 = [2 -2 2 0;
     -1 -6 0 -1;
     -4 -1 -1 -4;
     -5 0 -6 5];

A2 = [-3 7 -1 -1;
      -3 1 2 0;
      -1 -3 0 3;
       0 1 1 2];

A3 = [0 2 7 3;
      0 0 1 -1;
      0 0 -4 -1;
      -3 -1 2 0];
  
%%% LMI test
lmisolver = 'sedumi';
N = 100;
nu = logspace(-4,-1,N);

%% Using inv(A1)
V{1} = -inv(A1)*A2;
V{2} = -inv(A1)*A3;

for i = 1:N
   tsol(i) = remark3(V,nu(i),lmisolver);
end
figure,semilogx(nu,tsol)
ylim([-2 1])
title('Using inv(A1)')
xlabel('\nu')
ylabel('t_{min}')

%% Using inv(A2)
V{1} = -inv(A2)*A1;
V{2} = -inv(A2)*A3;

for i = 1:N
   tsol(i) = remark3(V,nu(i),lmisolver);
end
figure,semilogx(nu,tsol)
ylim([-2 1])
title('Using inv(A2)')
xlabel('\nu')
ylabel('t_{min}')

%% Using inv(A3)
V{1} = -inv(A3)*A1;
V{2} = -inv(A3)*A2;

for i = 1:N
   tsol(i) = remark3(V,nu(i),lmisolver);
end
figure,semilogx(nu,tsol)
ylim([-2 1])
title('Using inv(A3)')
xlabel('\nu')
ylabel('t_{min}')

disp('LMI test outcome:')
disp('If tmin = -1 for some \nu in any of the plots, A(alpha) will have full rank for all alpha in the unit simplex.')
disp('Otherwise, the test is inconclusive.')
disp('Press any key to continue...')
pause

%%% Cross-vertex test
eig(A1'*A2 + A2'*A1)
eig(A1'*A3 + A3'*A1)
eig(A2'*A3 + A3'*A2)

disp('Cross-vertex test outcome:')
disp('If all eigenvalues are positive, A(alpha) will have full rank for all alpha in the unit simplex.')
disp('Otherwise, the test is inconclusive.')
disp('Press any key to continue...')
pause

%%% P-matrix test
I = eye(4);
H = [A1*inv(A3) (A2 - A1)*inv(A3);
     -I I];
for i = 1:size(H,1)
    % Evaluating the principal minors of H
     det(H(1:i,1:i))
end

disp('P-matrix test outcome:')
disp('If all principal minors are positive, A(alpha) will have full rank for all alpha in the unit simplex.')
disp('Otherwise, the test is inconclusive.')

%%% Plot of det(A(alpha)) for different values of alpha in the unit simplex
rand('state',0)
Nreal = 10000;
for j = 1:Nreal
    aux1 = rand; aux2 = rand; aux3 = rand;
    sumaux = aux1 + aux2 + aux3;
    a1(j) = aux1/sumaux;
    a2(j) = aux2/sumaux;
    a3(j) = aux3/sumaux;
    A = a1(j)*A1 + a2(j)*A2 + a3(j)*A3;
    detA(j) = det(A);
end

% 3D plot with a1, a2 axes
plot3(a1,a2,detA,'k.')
box on
xlabel('\alpha_1')
ylabel('\alpha_2')
zlabel('det(A(\alpha))')
h = gca;set(h,'Xtick',[0 0.5 1],'Ytick',[0 0.5 1],'FontSize',12)
hold on, plot3([0 0],[0 0],[-250 0],'k')
plot3([0 0],[0 1],[0 0],'k')
plot3([0 1],[0 0],[0 0],'k')
view(-20,25), grid
patch([0 1 0],[0 0 1],[-250 -250 -250],[0.5 0.5 0.5])
zlim([-250 0])
