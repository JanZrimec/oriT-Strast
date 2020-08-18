function [out] = strast_inner_loopV6_2(target,querySet,idx,distMat,nq,tresh,fun,lnic,ival)
% 2_11_17: popravki in propagacija pval algo call
% 31_10_17
% remove outMin
% add all four orientations
% add pval treshold hit finding function

% popravek 28_2_16 - SEQRCOMPLEMENT
% add compensation to targetSet for 
% A. NNx processing - nepotreno ker enako za target in query
% B. length of orit edge to nic = cirkularizacija

lnic2 = length(querySet{1,1})-lnic;
target2 = [target(end-lnic+1:end),target,target(1:lnic2)];
targetRC = [target(end-lnic2+1:end),target,target(1:lnic)];
% targetR = targetRC; %ok
% targetC = target2;
out = []; 

for j = 1:nq
    
    % get alignment for orientation
    % get hits at treshold of pval
    F1 = onedalignV5_2(target2,querySet{j,1},idx,distMat);
    out{j,1} = getHitsV6(F1,tresh,fun);
    out{j,5} = getHitsV6_2(F1,tresh,fun,querySet{j,2}(1,ival));
    % RC
    F2 = onedalignV5_2(targetRC,seqcomplement(querySet{j,1}(end:-1:1)),idx,distMat);
    out{j,2} = getHitsV6(F2,tresh,fun);
    out{j,6} = getHitsV6_2(F2,tresh,fun,querySet{j,2}(1,ival));
    % R
    F3 = onedalignV5_2(targetRC,querySet{j,1}(end:-1:1),idx,distMat);
    out{j,3} = getHitsV6(F3,tresh,fun);
    out{j,7} = getHitsV6_2(F3,tresh,fun,querySet{j,2}(1,ival));
    % C
    F4 = onedalignV5_2(target2,seqcomplement(querySet{j,1}),idx,distMat);
    out{j,4} = getHitsV6(F4,tresh,fun);
    out{j,8} = getHitsV6_2(F4,tresh,fun,querySet{j,2}(1,ival));
%     F1o{j} = F1;
%     F1m{j} = mean(F1);
end

end