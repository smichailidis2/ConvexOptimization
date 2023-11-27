function [x,fun_val,f_traj,t_total,iter] = newton_method_backtracking(f,A,b,g,H, ...
    x0,s,alpha,beta,epsilon)
% Newton method with backtracking stepsize rule
%
%=======================================
% ------ INPUTS ------
% f ......... objective function
% g ......... gradient of the objective function
% H ......... hessian of the objective function
% x0......... initial point
% s ......... initial choice of stepsize
% alpha ..... tolerance parameter for the stepsize selection
% beta ...... the constant in which the stepsize is multiplied
% at each backtracking step (0<beta<1)
% epsilon ... tolerance parameter for stopping rule
%
% A,b ......  parameters of logarithmic function. Set to 'none' if 
%
% ------ OUTPUTS ------
% x ......... optimal solution (up to a tolerance)
% of min f(x)
% fun_val ... optimal function value
% iter ........ total iterations
%=======================================

x=x0;
grad=g(x);
fun_val=f(x);
hessian=H(x);
iter=0;

f_traj =[];
f_traj(end+1,:) = fun_val;

lambda = 0.5*grad'*hessian*grad;

tic;

while (lambda>epsilon)
    iter=iter+1;
    t=s;
    %newton step
    DX = hessian\grad;

    % feasability check
    while (~domain_check(A,b,x-t*DX))
        t=beta*t;
    end

    while (fun_val-f(x-t*DX)<alpha*t*grad'*DX)
        t=beta*t;
    end

    x = x - t*DX;

    fun_val=f(x);
    grad=g(x);
    hessian=H(x);
    %fprintf('iter_number = %3d norm_grad = %2.6f fun_val = %2.6f \n',iter,norm(grad),fun_val);
    lambda = 0.5*grad'*hessian*grad;

    f_traj(end+1,:) = fun_val;
end
t_total = toc;


end