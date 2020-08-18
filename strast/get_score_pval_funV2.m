function x=get_score_pval_funV2(y,p)
% 6_11_2017: fit only to point 10e-1 (4th parameter in p1), above is const.
% thus assert that x < p(4)

x = (-p(2)+sqrt(p(2).^2-4.*p(1).*(p(3)-y)))./(2.*p(1)); % this is always on the
% left side of the function so x < p(4)
x = real(x);
% if x > p(4) % failsafe 
%     x = p(4);
% end
x(x > p(4)) = p(4); % vectorized

end