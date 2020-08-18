function [out,err] = strast_mainV6_2(targetSet,querySet,idxMAT,dMAT,vals,pMAT,vals2,tresh,lnic)
%% 9.11.17 
% - optimization to enable finding all possible new orit in data
% - optimization of certain components
% - optimization of output files for saving

% input:
% - target, query sets
% - idx and distmat whole matrices
% - values at which idx and distmat parameters to search (multuple)
% - lnic - nic site location in sequence after n bp from 5' to 3' (offset)

% output:
% - out: 
% (4 alignment vectors per input sequence, per each of four possible 
% oritentations of input: cell of {target,query}{4xarrays of
% uint32=4bytes}) == too much data
% - matrices of hits below treshold values -(6,7,8,9,10,12,14,16,18,20)
% - outMin - remove searching for minimum with other output processing algo

% ACGT or nonACGT:
% - combine both algorithms 
% - expand nonACGT algorithm with ACGT to speed up - locate nonACGT points
% - include ACGT determination directly in algorithm - no input needed

% data structure and files:

% removal:
% - (pvalues calculation and tresholds) NE

% options:
%  - input p-val functions and tresholds to find hits below tresholds (-6,
%  -10) - JA
% - bag all hits - JA
% - propp - pvalues for propagation of putative nic according to
% experimental, at the different values of parameters given

%% cycle thorugh parameters
% cycle through targetSet and through querySet
nt = size(targetSet,1);
nq = size(querySet,1);
err=[];
out = cell(nt,size(vals,1)); % cell output je v redu glede prostora - pomembne so 
% interne spremenljivke, zelo hitro se isce cez njega v primerjavi struct
% inside each (nt,nq) is (9,4) cell with list of hits for 4 orientation at 
% 9 treshold values

for i=1:size(vals,1)
    idx = idxMAT{vals(i,1),vals(i,2)};
    distMat = dMAT{vals(i,1),vals(i,2)};
    fun = pMAT{vals2(i,1),vals2(i,2)};
    lnum=lnic(i);
    
    parfor ii = 1:nt
        try
            disp(ii)
        % if isNotATCG(targetSet{ii,1})
        %    [out{ii,i}] = strast_inner_loopV5_1(targetSet{ii,1},querySet,idx,distMat,nq,tresh,fun,lnic);     
        % else
        %    [out{ii,i}] = strast_inner_loopV5(targetSet{ii,1},querySet,idx,distMat,nq,tresh,fun,lnic);
        % end
        [out{ii,i}] = strast_inner_loopV6_2(targetSet{ii,1},querySet,idx,distMat,nq,tresh,fun,lnum,i);  
        catch
            err=[err,ii];
        end
    end
end

end