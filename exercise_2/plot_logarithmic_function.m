function plot_logarithmic_function(f,A,b,x_opt)
%
% f  -----> function handle
% A,b ----> parameters of f
% x_opt --> optimal point of f
%
samples = 128;
m = length(b);

x_1 = linspace(x_opt(1) - 1, x_opt(1) + 1, samples);
x_2 = linspace(x_opt(2) - 1, x_opt(2) + 1, samples);

l1 = length(x_1);
l2 = length(x_2);

f_x1_x2 = zeros(l1,l2);

for i = 1:l1
    for j = 1:l2
        x_ij = [x_1(i) ; x_2(j)];
        
        % If the element does not belong in the domain of f value of
        % f(x1,x2) is set to inf.
        if ( sum( (b - A*x_ij) < zeros(1,m) ) > 0 )
            f_x1_x2(i,j) = inf;
        else
            f_x1_x2(i,j) = f(x_ij);
        end
    end
end

%details(f_x1_x2)
%class(f_x1_x2)

[X1,X2] = meshgrid(x_1,x_2);

% mesh
figure('Name','Logarithmic function near optimal x')
sc = meshc(X1,X2,f_x1_x2,'FaceColor','flat');
colormap("winter");
colorbar;
sc(2).EdgeColor = '#00A36C';
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
zlabel('$f(\mathbf{x}) = \sum_{i=1}^{m}log(b_i - a_i^T\mathbf{x}) $','Interpreter','latex')
grid on;
hold on;
plot3(x_opt(1),x_opt(2),f(x_opt),'bo','MarkerSize',5,'MarkerFaceColor',[1, 0, 0])
legend('$f(x_1 , x_2)$','level sets','$f(\vec x_*)$','Location','best','Interpreter','latex')

% contours 
figure('Name','Contour plot of function near optimal x')
contour(X1,X2,f_x1_x2)
colormap("summer");
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
grid on;
hold on;
plot3(x_opt(1),x_opt(2),f(x_opt),'bo','MarkerSize',5,'MarkerFaceColor',[1, 0, 0])
title('Level sets of $f$','Interpreter','latex')
legend('','$f(\vec x_*)$','Location','best','Interpreter','latex')

fprintf("\nOptimal point shown in red\n\n")


end