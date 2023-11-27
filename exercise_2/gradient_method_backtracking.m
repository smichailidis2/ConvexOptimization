function [x,fun_val,traj,f_traj,t_total,iter]=gradient_method_backtracking(f,g,x0,s,alpha,beta,epsilon,func_type,A,b)
% Gradient method with backtracking stepsize rule
%
%=======================================
% ------ INPUTS ------
% f ......... objective function
% g ......... gradient of the objective function
% x0......... initial point
% s ......... initial choice of stepsize
% alpha ..... tolerance parameter for the stepsize selection
% beta ...... the constant in which the stepsize is multiplied
% at each backtracking step (0<beta<1)
% epsilon ... tolerance parameter for stopping rule
%
% func_type .... logarithmic or quadratic (string)
% A,b ......  parameters of logarithmic function. Set to 'none' if 
%
% ------ OUTPUTS ------
% x ......... optimal solution (up to a tolerance)
% of min f(x)
% fun_val ... optimal function value
%
% (my additions)
% traj ....... trajectory of point x
% t_total ..... total time of execution
% iter ........ total iterations
%=======================================

x=x0;
grad=g(x);
fun_val=f(x);
iter=0;

traj = [];
f_traj = [];
traj(end+1,:) = x;
f_traj(end+1,:) = fun_val;
tic;

if ( strcmp( func_type,'quadratic') )
    while (norm(grad)>epsilon)
        iter=iter+1;
        t=s;
        while (fun_val-f(x-t*grad)<alpha*t*norm(grad)^2)
            t=beta*t;
        end
        x=x-t*grad;
        fun_val=f(x);
        grad=g(x);
        %fprintf('iter_number = %3d norm_grad = %2.6f fun_val = %2.6f \n',iter,norm(grad),fun_val);
        traj(end+1,:) = x;
        f_traj(end+1,:) = fun_val;
    end
elseif ( strcmp(func_type,'logarithmic') )

    if( strcmp(A,'none') && strcmp(b,'none')  )
        error("\n\tgradient_method_backtracking.m: Incorrent input arguments <A>/<b>")
    end
    while (norm(grad)>epsilon)
        iter=iter+1;
        t=s;

        % feasability check
        while (~domain_check(A,b,x-t*grad))
            t=beta*t;
        end
        %fprintf("t after feasibility check = %3d\n",t)

        while (fun_val-f(x-t*grad)<alpha*t*norm(grad)^2)
            t=beta*t;
        end
        x=x-t*grad;

        fun_val=f(x);
        grad=g(x);
        %fprintf('iter_number = %3d norm_grad = %2.6f fun_val = %2.6f \n',iter,norm(grad),fun_val);
        traj(end+1,:) = x;
        f_traj(end+1,:) = fun_val;
    end
else
    error("\n\tgradient_method_backtracking.m: Incorrent input argument <func_type>");
end

%fprintf("\nOptimal function value = %d\n", f(x))

t_total = toc;

end

