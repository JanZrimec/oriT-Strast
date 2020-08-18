function out=getHitsV6_2(in,tresh,fun,propp)
% 2_11_17: upgrade to include propagation of pvals in finding hits
% 2 option: 1. score to pval 2. pval to score
% 31_10_17
% list of hits below treshold values additive from lowest to highest
% disregard around best result in window of +/-wind bp

wind = 5;
out = cell(length(tresh),1);
% mf = max(in);

% transform score vector to pval vector

% propagate all pvals 
% p(A or B) = p(A) + p(B) - p(A and B)
% p(A and B) = p(A) * p(B)

pval = get_pval_score_funV2(in,fun);
pval = log10(10.^pval+10.^propp-10.^(pval+propp));
mf = max(pval);

for i=length(tresh):-1:1 % filamo od zadaj brez ponavljanj hitov

    out{i,1} = single(find(pval < tresh(i))); % razlicno stevilo elementov
    
    if ~isempty(out{i,1})
        for k=1:length(out{i,1}) % za vsakega shranim tudi pravi pval
            out{i,1}(k,2) = pval((out{i,1}(k,1)));
        end
        
        % delete window - pazi robove
        aaa = (out{i,1}(:,1));
        for k=1:wind
            aaa(:,k+1) = aaa(:,1)+k;
            aaa(:,2*wind+2-k)=aaa(:,1)-k;
        end
        aaa2 = unique(aaa(:));
        
    pval(aaa2((aaa2)>0&(aaa2)<=length(pval))) = mf; % delete %in(out{i,1}(:,1)) = mf;
    end
end

end

% aaa=[100,1000,10000];
% for i=1:wind
%     aaa(i+1,:)=aaa(1,:)+i;
%     aaa(2*wind+2-i,:)=aaa(1,:)-i;
% end
% aaa(:)