function [x,fun_val,traj,f_traj,t_total,iter] = gradient_exact_method(f,g,P,x0,epsilon)
% Gradient method with exact stepsize rule
%
%=======================================
% ------ INPUTS ------
% f ......... objective function
% g ......... gradient of the objective function
% x0 ......... initial point
% P ......... p.s.d matrix of 
% % quadratic input function f
% epsilon ... tolerance parameter for stopping rule
%
% ------ OUTPUTS ------
% x ......... optimal solution (up to a tolerance)
% of min f(x)
% fun_val ... optimal function value
%
% (my additions)
% traj ....... trajectory of point x
% f_traj ...... values of f for 
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
while (norm(grad)>epsilon)
    iter=iter+1;
    % t = argmin f( x - t*grad[f(x)] ) , t>=0
    t = ( norm(grad)^2 )/( grad'*P*grad );
    x=x-t*grad;
    fun_val = f(x);
    grad = g(x);
    %fprintf('iter_number = %3d norm_grad = %2.6f fun_val = %2.6f \n',iter,norm(grad),fun_val);
    traj(end+1,:) = x;
    f_traj(end+1,:) = fun_val;
end

t_total = toc;

end


