function H = log_hessian(A,b,x)
% Hessian of logarithmic function

H = 0;

for i=1:length(b)
    H = H + ( A(i,:)'*A(i,:) )/( b(i) - A(i,:)*x )^2;
end

end