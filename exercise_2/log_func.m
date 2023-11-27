function f = log_func(x,c,b,A)

        f = c'*x - sum(log(b-A*x));
        
end