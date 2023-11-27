%%
% EXERSISE 1 - OPTIMIZATION
%
% MICHAILIDIS STERGIOS 2020030080
%
% winter 2023
%%
close all
clear
clc
% for better text, change text interpreter to latex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
%%
% 1)

% x axis
x = linspace(0,5, 20);
% taylor point
x_0 = 3;

% function
f = 1./(x + 1);

figure(1)
plot(x,f)
hold on;

% first order taylor approx.
f1 = 1/(1 + x_0) - (x - x_0)./(1 + x_0)^2 ;
plot(x,f1)

% second order taylor approx.
f2 = f1 + 0.5*( (x - x_0).^2 )./(1 + x_0)^3 ;
plot(x,f2)

hold off;
grid on;
xlabel('x axis')
ylabel('y axis')
title('$f(x) , f_{(1)}(x), f_{(2)}(x)$')
legend('f(x)', 'First order Taylor', 'Second order Taylor')

%%
% 2)
% a.  d.  e.

% axes
x_1 = x;
x_2 = x;
len = length(x);

% taylor point
x_0 = [3,3];

% initialize the dunction and its first and second taylor approx.
f = zeros(len, len);
f1 = f;
f2 = f1;

for j = 1:len
    for k = 1:len
        % function 
        f(j,k) = 1/(1 + x_1(j) + x_2(k));
        
        % first order taylor approx.
        f1(j,k) = 1/(1 + x_0(1) + x_0(2)) - ([x_1(j) - x_0(1) , x_2(k) - x_0(2)])*[1/(1 + x_0(1) + x_0(2))^2; 1/(1 + x_0(1) + x_0(2))^2];
        
        % second order taylor approx.
        f2(j,k) = f1(j,k) + 0.5*([x_1(j) - x_0(1) , x_2(k) - x_0(2)]) * (1 + x_0(1) + x_0(2))^(-3) * [2 2; 2 2] * ([x_1(j) - x_0(1); x_2(k) - x_0(2)]);
            
        
    end
end

figure(2)
mesh(x_1, x_2, f,FaceColor="texturemap");
hold on;
mesh(x_1,x_2, f1,EdgeColor = [0.8 0.5 0.2],FaceColor = [0.8 0.5 0.2])
mesh(x_1,x_2, f2,EdgeColor = [0 0 0],FaceColor = [0.24 1 0])

xlabel('$x_1$')
ylabel('$x_2$')
axis tight
zlabel('$f(x_1,x_2)$ , first and second order approx.')
legend('f(x)', 'First order Taylor', 'Second order Taylor')
hold off;

% b.

figure(3)
contourf(x_1, x_2, f)
grid on;
xlabel('$x_1$')
ylabel('$x_2$')
title('Contour plot of $$ f(x_1,x_2) = \frac{1} {1+x_1 + x_2} $$');
%%
% 5)
clear i
% c.
a = -3 : 1 : 3;
f = zeros(length(x), length(a));
figure(4)
for i = 1 : length(a)
    f(:,i) = x.^a(i);
    plot(x, f(:,i), 'LineWidth', 1.5)
    xlabel('x axis')
    ylabel('f(x)')
    hold on
end
grid on
title('$f(x) = x^{a}$ for different values of $\textit{a}$')
legend('$\textit{a} = -3$','$\textit{a} = -2$', '$\textit{a} = -1$', '$\textit{a} = 0$', '$\textit{a} = 1$', '$\textit{a} = 2$', '$\textit{a} = 3$')
axis([x(1) x(end) x(1) x(end)])


% d.

x = -10:0.2:10;
% Euclidean norm and Euclidean norm^2
x_1 = x;
x_2 = x;
len = length(x);

f = zeros(len, len);
for i = 1 : len
    for j = 1 : len
        f(i,j) = x_1(i)^2 + x_2(j)^2;
    end
end

f_n = sqrt(f);

% plots
figure(5)
mesh(x_1,x_2,f_n)
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\|\mathbf{\vec{x}}\|_2$')
title('Euclidian norm')
grid on;


figure(6)
mesh(x_1,x_2,f)
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$\|\mathbf{\vec{x}}\|_2^2$')
title('Euclidian norm squared')
grid on;

%%
clear
% generate radom matrix P, vector q and scalar r:
A = randn(2,2);
P = A*A';
q = randn(2,1);
r = randn;

x_star = -inv(P)*q
f_min = 0.5*x_star'*P*x_star + q'*x_star + r

x1 = x_star(1)-10:0.2:x_star(1)+10;
x2 = x_star(2)-10:0.2:x_star(2)+10;

[x_1,x_2] = meshgrid(x1, x2);

for i = 1:size(x_1)
    for j = 1:size(x_1)
        x_t = [x_1(i,j); x_2(i,j)];
        f(i,j) = 0.5*x_t'*P*x_t + q'*x_t + r;
    end
end

% mesh plot
figure(7)
mesh(x_1,x_2,f);
grid on;
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$f(\vec{x})$')
title('Plot $f(\vec{x}) = \frac{1}{2}\vec{x}^T\mathbf{P}\vec{x} + q^T \vec{x} + r$')
hold on;
plot(x_star(1), x_star(2), '.r')

% contour plot
figure(8)
contourf(x_1,x_2,f);
grid on;
xlabel('$x_1$')
ylabel('$x_2$')
title('Contour of $f(\vec{x}) = \frac{1}{2}\vec{x}^T\mathbf{P}\vec{x} + q^T \vec{x} + r$')
hold on;
plot(x_star(1), x_star(2), '.r')

