function g = log_gradient(x,c,b,A)
% 
% outputs the gradient of logarithmic function  
%

g = c;
for i=1:length(b)
    g = g + A(i,:)'/(b(i) - A(i,:)*x);
end

end