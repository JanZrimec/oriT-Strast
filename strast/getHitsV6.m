function out=getHitsV6(in,tresh,fun)
% 31_10_17
% list of hits below treshold values additive from lowest to highest
% disregard around best result in window of +/-wind bp

wind = 5;
out = cell(length(tresh),1);
mf = max(in);

for i=length(tresh):-1:1 % filamo od zadaj brez ponavljanj hitov

    score = get_score_pval_funV2(tresh(i),fun);
    out{i,1} = single(find(in < score)); % razlicno stevilo elementov
    
    if ~isempty(out{i,1})
        for k=1:length(out{i,1}) % za vsakega shranim tudi pravi pval
            out{i,1}(k,2) = get_pval_score_funV2(in(out{i,1}(k,1)),fun);
        end
        
        % delete window - pazi robove
        aaa = (out{i,1}(:,1));
        for k=1:wind
            aaa(:,k+1) = aaa(:,1)+k;
            aaa(:,2*wind+2-k)=aaa(:,1)-k;
        end
        aaa2 = unique(aaa(:));
        
    in(aaa2((aaa2)>0&(aaa2)<=length(in))) = mf; % delete %in(out{i,1}(:,1)) = mf;
    end
end

end

% aaa=[100,1000,10000];
% for i=1:wind
%     aaa(i+1,:)=aaa(1,:)+i;
%     aaa(2*wind+2-i,:)=aaa(1,:)-i;
% end
% aaa(:)