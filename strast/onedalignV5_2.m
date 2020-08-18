function [F] = onedalignV5_2(target,query,idx,distMat)
% 1_11_17: add if to have one function
% 1_11_17: expand to optimized version of onedalignV3 for nonACGT
% 1_11_17: vectorized version of onedalignV2_2_speed
% 10_4_17: upgrade of onedaling 2_9_16 using precalculated distance matrix

N = log10(length(idx))/log10(4);
nt = length(target)-N+1; 
nq = length(query)-N+1;
F = zeros(nt-nq+1,1);
[intseq1,intout2] = nt2intV2_2(target); % position of odd chars are stored
intseq2 = nt2intV2_1(query);

% create vectorization matrix 1
mat = [];
for i = 1:N
    mat = [mat, (i:nt-1+i)'];
end
    
% get data
numseq1 = base2dec(intseq1(mat)-1)+1;
numseq2 = base2dec(intseq2(mat(1:nq,:))-1)+1;

% create blacklist - all NNxs that include odd chars
if ~isempty(intout2)
    tmp = length(intout2);
    intout = zeros(N+N-1,tmp);
    for i = 1:tmp
        intout(:,i) = (intout2(i)-N+1:intout2(i)+N-1)';
    end

% at odd positions distance is max - expand matrices to include max dist
    numseq1(intout) = length(idx)+1;
    idx(end+1) = length(distMat)+1;
    distMat(end+1,:) = max(max(distMat));
end

% get F
for i=1:nq % make for in smallest dim nq
    % sum diagonally over all target points in nt vs nq
    F = F + distMat(idx(numseq1(i:nt-nq+i)),idx(numseq2(i)));
end

end