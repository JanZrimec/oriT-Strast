function x=funs_polynomial_rev(y,p)

%y = p(1).*x.^2 + p(2).*x + p(3);

x(1) = (-p(2)+sqrt(p(2)^2-4*p(1)*(p(3)-y)))/(2*p(1)); %y gre v p3
%x(2) = (-p(2)-sqrt(p(2)^2-4*p(1)*(p(3)-y)))/(2*p(1));

%x = solve(p(1)*x^2 + p(2)*x + p(3) - y == 0,x); 