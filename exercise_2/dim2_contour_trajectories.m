function dim2_contour_trajectories(P,q,traj, iter, str_algorithm_used)
% dim2_contour_trajectories
% ===========================
% INPUTS:
% P ........ p.s.d. matrix of quadratic
% function f = 0.5x' * P * x + q' * x
% q ........ 2x1 real vector of f
% traj ....... trajectory of minimal value x_*
% iter ....... number of iterations
% str_algorithm_used .... type of
% algorithm used (Backtracking, exact line search e.t.c.)
% string format
% ===========================


samples = 100;

x1 = linspace(min(traj(:,1)) - 1, max(traj(:,1)) + 1 , samples);
x2 = linspace(min(traj(:,2)) - 1, max(traj(:,2)) + 1 , samples);

% find the exact contours of the trjaectory of x_*
contours = zeros(1,iter);
for k = 1:iter
    contours(1,k) = 0.5*traj(k,:)*P*traj(k,:)' + q'*traj(k,:)';
end

[X1,X2] = meshgrid(x1,x2);

for i = 1:length(X1)
    for j = 1:length(X2)
        f(i,j) = 0.5*[X1(i,j);X2(i,j)]'*P*[X1(i,j);X2(i,j)] + q'*[X1(i,j);X2(i,j)];
    end
end

F = reshape(f,size(X1));

txt = str_algorithm_used + " with " + iter + " iterations";
figure('Name',str_algorithm_used)
if (iter > 1) 
    contour(X1,X2,F,contours)
elseif (iter == 1)
    contour(X1,X2,F)
end
hold on;
plot(traj(:,1),traj(:,2),'-ro','MarkerFaceColor','r',LineWidth=0.9)
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title(txt);
legend('level sets','trajectory of $x_*$','Location','northeastoutside')
grid on;
axis tight
hold off;


end