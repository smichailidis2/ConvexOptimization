%% OPTIMIZATION EXERCISE 2 %%
% 
% MICHAILIDIS STERGIOS 
% 2020030080
%
% November 2023
%
%%
clear
close all;
clc

set (groot ,'defaulttextinterpreter','latex');
set (groot ,'defaultAxesTickLabelInterpreter','latex');
set (groot ,'defaultLegendInterpreter','latex');

%format("long");

%% ====== Question B. ======
% (a)

% dimention number ; interchangable ; n>=2
n = 15;

A = randn(n,n);

[U,S,V] = svd(A);

observe_1 = U*U'
observe_2 = U'*U
fprintf('Observation: Matrix U is unitary!\nPress any key to continue.') ,pause;


% (b)
%% condition number 
K = 10;
%%
% construct the eigenvalues
l_min = rand();
l_max = K*l_min;

z = l_min + (l_max - l_min)*rand((n-2),1);

eig_P = [l_min;(z);l_max];

L = diag(eig_P);

% spectral decomposition of random matrix P
P = U*L*U';

% create random vector q
q = randn(n,1);

%% -- CLOSED FORM SOLUTION --

% x_star = argmin f(x) <=> nabla [f(x_star)] = 0 <=> x_star = -P^(-1) * q
x_star = -P\q

% p_star = f(x_star)
p_star = (1/2)*x_star'*P*x_star + q'*x_star


% --CVX--
cvx_begin
    
    variables x(n)
    minimize ( (1/2)*x'*P*x + q'*x )

cvx_end

fprintf('Optimal value (closed form) = %.5f\n',p_star)
fprintf('\nWe can see that the closed form solution is basically the same as the cvx solution...\n')

%% - GRADIENT DESCENT -
f = @(x) (1/2)*x'*P*x + q'*x;
g = @(x) P*x + q;

% required parameters
x0 = 2*randn(n,1);
s = 2;
alpha = 0.1;
beta  = 0.7;
epsilon = 10^-5;

% check the parameters
str = "Backtracking algorithm";
check_integrity(s,alpha,beta,epsilon,str)

% first lets use the backtracking method
[x,fun_val,traj,f_traj,t_total,iter]=gradient_method_backtracking(f,g,x0,s,alpha,beta, ...
    epsilon,'quadratic','none','none');

fprintf('\nBacktracking algorithm completed\n')
fprintf('\n\tTotal time = %f \n\tTotal number of iterations = %i\n\n',t_total,iter)
fprintf('\noptimal value of x* is :\n')
disp(x)
fprintf('\noptimal value of p* is :\n')
disp(f(x))

if(n<3)
    dim2_contour_trajectories(P,q,traj, iter, str)
end

% secondly we have th exact line search method
epsilon2 = 10^-5;

% check the parameters
str2 = "Exact line algorithm";
check_integrity('none','none','none',epsilon2,str2)

[x2,fun_val2,traj2,f_traj2,t_total2,iter2] = gradient_exact_method(f,g,P,x0,epsilon2);

fprintf('\nExact line algorithm completed\n')
fprintf('\n\tTotal time = %f \n\tTotal number of iterations = %i\n\n',t_total2,iter2)
fprintf('\n optimal value of x* is :\n')
disp(x2)
fprintf('\n optimal value  p* is :\n')
disp(f(x2))

if(n<3)
    dim2_contour_trajectories(P,q,traj2, iter2, str2)
end

%% log error plots

% initialize axis of iterations k
k_1 = 0:1:iter;
k_2 = 0:1:iter2;
% compute the log error for both gradient algorithms
L_k_backtrack = log(f_traj - p_star);
L_k_exact = log(f_traj2 - p_star);

% plot
figure('Name','Log-error plot')
plot(k_1,L_k_backtrack,'LineWidth',2)
hold on;
plot(k_2,L_k_exact,'LineWidth',1.5)
grid on;
legend('backtracking','exact')
xlabel('iterations $\mathbf k$')
ylabel('$log(f(x_k) - p_*)$')
title('Logarithmic error value vs. number of iterations')


%% - Convergence analysis -
% turn off for speed
flag = 1;
if(flag)
    fprintf('\nInitiate Convergence analysis. press any key to continue.\n'), pause;
    
    loops = 1000;
    % error tolerance parameter
    err = 10^-4;
    % dimention
    n = 2;
    % initialize hyperparameters
    alpha = 0.2;
    beta = 0.7;
    s = 2;
    
    % iteration numbers
    % actual
    iter_act_bt = 0;
    iter_act_exact = 0;
    % estimated
    iter_est_bt = 0;
    iter_est_exact = 0;
    
    str = "exercise2.m: Convergence analysis";
    
    for i=1:loops
    
        A = randn(n,n);
        [U,~,~] = svd(A);
    
        l_min = rand();
        l_max = K*l_min;
        
        z = l_min + (l_max - l_min)*rand((n-2),1);
        
        eig_P = [l_min;sort(z);l_max];
        
        L = diag(eig_P);
        
        % spectral decomposition of random matrix P
        P = U*L*U';
        
        % create random vector q
        q = randn(n,1);
    
        x0 = randn(n,1);
        
        f = @(x) (1/2)*x'*P*x + q'*x;
        g = @(x) P*x + q;
    
        % closed form solutions
        x_star = -inv(P)*q;
        p_star = (1/2)*x_star'*P*x_star + q'*x_star;
        
        check_integrity(s,alpha,beta,err,str)
        % gradient methods
        [~,fun_val_bt,~,~,~,iter_bt]=gradient_method_backtracking(f,g,x0,s,alpha,beta, ...
        err,'quadratic','none','none');
        [~,fun_val_exact,~,~,~,iter_exact] = gradient_exact_method(f,g,P,x0,err);
        
        iter_act_bt = iter_act_bt + iter_bt;
        iter_act_exact = iter_act_exact + iter_exact;
    
        % |estimations for current loop
        % exact
        est_exact = abs (ceil( K*abs( log( (fun_val_exact - p_star)/err ) ) ));
        iter_est_exact = iter_est_exact + est_exact;
        % backtracking
        c = 1 - min( 2*l_min*alpha, (2*beta*alpha*l_min)/l_max );
        est_bt = abs( ceil( log( (fun_val_bt - p_star)/err )/log(1/c) ) );
        iter_est_bt = iter_est_bt + est_bt;
    end
    
    average_iter_bt_estimated = iter_est_bt/loops;
    average_iter_exact_estimated = iter_est_exact/loops;
    
    average_iter_bt_actual = iter_act_bt/loops;
    average_iter_exact_actual = iter_act_exact/loops;
    
    fprintf("\n=======================================================================\n")
    fprintf("\n\tAverage backtracking iterations = %.2f\n",average_iter_bt_actual)
    fprintf("\n\tAverage exact line iterations = %.2f\n",average_iter_exact_actual)
    fprintf("\n|Estimations:|\n")
    fprintf("\n\tAverage number of backtracking ietrations (estimated) = %.2f\n",average_iter_bt_estimated)
    fprintf("\n\tAverage number of exact line ietrations (estimated) = %.2f\n",average_iter_exact_estimated)
    fprintf("\n=======================================================================\n")

end
%%
fprintf('\nMoving on to question C. press any key to continue.\n'), pause;
fprintf('\nClearing workspace variables...')
clearvars;

%% ====== QUESTION C. ======

% matrix and vector dimentions
n = 20;
m = 300;

c = rand(n,1);
b = rand(m,1);
A = randn(m,n);

f = @ (x) log_func(x,c,b,A);
g = @ (x) log_gradient(x,c,b,A);

%% (a)
% -- CVX --
cvx_begin
    
    variables x(n,1)
    minimize ( (c')*x - sum( log(b - A*x) ) )

cvx_end

%% (b)
% plot the logarithmic
TF = isinf(cvx_optval);
if ( n == 2 && ~TF)
    plot_logarithmic_function(f,A,b,x);
end

%% (c)
% backtracking descent for logarithmic
x0 = zeros(n,1);
%t_0 = 1
s = 1;
alpha = 0.1;
beta  = 0.4;
epsilon = 10^-3;
% check the parameters
str = "Backtracking algorithm";
check_integrity(s,alpha,beta,epsilon,str)

if (~TF)
[x_gm,fun_val,traj,f_traj,t_total,iter]=gradient_method_backtracking(f,g,x0,s,alpha,beta, ...
    epsilon,'logarithmic',A,b);

fprintf("\n=========================================================\n")
fprintf('\nBacktracking algorithm completed\n')
fprintf('\n\tTotal time = %f \n\tTotal number of iterations = %i\n\n',t_total,iter)
fprintf('\nAnd optimal value of x* is :\n')
disp(x_gm)
fprintf("\nOptimal function value = %.4f\n", f(x_gm))
end

%% (d) - Newton method
% Declare the hessian
H = @ (x) log_hessian(A,b,x);

% check the parameters
str = "Backtracking algorithm";
check_integrity(s,alpha,beta,epsilon,str)

if (~TF)
    [x_nt,func_val_nt,f_traj_nt,t_total_nt,iter_nt] = newton_method_backtracking(f,A,b,g,H,x0, ...
        s,alpha,beta,epsilon);

fprintf("\n=========================================================\n")
fprintf('\nNewton Backtracking algorithm completed\n')
fprintf('\n\tTotal time = %f \n\tTotal number of iterations = %i\n\n',t_total_nt,iter_nt)
fprintf('\nAnd optimal value of x* is :\n')
disp(x_nt)
fprintf("\nOptimal function value = %.4f\n", f(x_nt))
end

%% (e) semilog error plots

if (~TF)
    % initialize axis of iterations k
    k_1 = 0:1:iter;
    k_2 = 0:1:iter_nt;
    
    % plot
    figure('Name','semilog error plot')
    semilogy(k_1,( f_traj-cvx_optval ),'LineWidth',2)
    hold on;
    semilogy(k_2,( f_traj_nt-cvx_optval ),'LineWidth',1.5)
    grid on;
    legend('Gradient backtracking','Newton backtracking')
    xlabel('iterations $\mathbf k$')
    ylabel('$log(f(x_k) - p_*)$')
    title('error value vs. number of iterations (log scale y axis)')
end
