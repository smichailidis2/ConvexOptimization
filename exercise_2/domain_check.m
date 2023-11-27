function flag = domain_check(A,b,x)

flag = 1;
for i = 1:length(b)
    if( b(i) - A(i,:)*x <= 0 )
        flag = 0;
        break;
    end
end

end