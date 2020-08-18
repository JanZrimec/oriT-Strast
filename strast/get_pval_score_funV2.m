function [out1] = get_pval_score_funV2(score,p1)
% 6_11_2017: fit only to point 10e-1 (4th parameter in p1), above is const.

% if score < p1(4)
%     out1 = funs_polynomial(score,p1);
% else
%     out1 = -1;
% end
% ne deluje za vektorje, for loop je pocasen
% vectorized
out1 = funs_polynomial(score,p1);
out1(score >= p1(4)) = -1; %~(score < p1(4))

end