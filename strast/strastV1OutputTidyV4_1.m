function [out,outcell,out2,outcell2] = strastV1OutputTidyV4_1(in1,cutoff,qtable) 
            % ,outcell2,hitcnt
% processes strast output to csv ready output of plasmid-orit entries
% results in much smaller file size

% 6_6_19  fix cutoff, orientation problems
% 23_7_18 fix filtering of unique locations within tolerance       

% 22_7_18 return {Loc, orient, query, mob, submob, pval}, -corrected pval
%         return everything found at all locations - unfiltered complete
%         return also only filtered best found
%         fix output to discernible table format   

% 20_6_18 fix obtaining from lower levels tha cutoff,
%         fix/add removal of 'duplicate' variants (albeit diff orient.)
%         (optional later) filter hits with diff relaxases

% in1 ..    single out file of strast - has uncorrected pvals (col 1-4)
% in2 -> removed
% out ..    each row is target with entries per location found
%           each entry is cell of {loc,orientation,query,pval}
% outcell.. returns tidy data that can be written to csv
% out2 ..   filtered out to keep in a certain window only best hit per MOB 
% outcell2..writeable
% (hitcnt .. how many hits per entry per target in a MOB group)-> TODO
% cutoff .. p-value at which to look (see tresh as input to strast)
%           tresh=-[6,7,8,9,10,12,14,16,18,20,30];
% querytable 
%        .. load file load('dataset_nic_112full_1_11_17.mat', 'NIC112')

% comment: from this a data structure is constructed and saved, or
% separate outputs are added to file

n = length(in1);
orient={'F','RC','R','C'};
out = cell(n,1);
outcell = cell(n,1);
out2 = cell(n,1);
outcell2 = cell(n,1);
%hitcnt = cell(n,1);
%tolerance = 1; % pval corrected tolerance
tolerance2 = 30; % same relaxase hit tolerance
cutoff1 = 1; % transform to scale used in script

for i=1:n % per plasmid
    
    if ~isempty(in1{i}) % errors in alginment run
    for qq=1:112 % per query
        for j=1:4 % per orientation
            for cc = length(cutoff):-1:cutoff1 % per cutoff level
             if ~isempty(in1{i}{qq,j}{cc})
                for k=1:length(in1{i}{qq,j}{cc}(1,1)) % multiple hits
                    
                    % make cell array with stored data
                    % test cell {2,1}{11,5}{1}(1:2)
                    tmpdata{1} = in1{i}{qq,j}{cc}(1,1); % 
                    tmpdata{2} = orient{j};
                    tmpdata{3} = qq;
                    tmpdata{4} = qtable{qq,2}; % MOB
                    tmpdata{5} = qtable{qq,3}; % subgroup
                    tmpdata{6} = in1{i}{qq,j}{cc}(1,2);
                    
                    out{i}(end+1,:) = tmpdata;
                end
             end
            
            end
            
        end
    end
    %out{i}(1,1) = []; % delete empty cell?
    
    %% filtering
    % check hits that are same element variants within tolerance2
    % in a MOB group, if different MOB group treat as separate
    % FIX: diff algorithm
    
    % Step 1: precalculate upper triangular matrix of comparisons and mark 
    % below tolerance
    tolerance_matrix = cell(size(out{i},1));
    for ii=1:size(out{i},1)-1
        for iii=ii+1:size(out{i},1)
            tolerance_matrix{ii,iii} = ...
                abs(out{i}{ii,1}-out{i}{iii,1}) <= tolerance2;
        end
    end
    % Step 2: move around upper triangle of matrix, when true check p-val,
    % remove higher from list of indexes
    list_removed = [];
    for ii=1:size(out{i},1)-1
        for iii=ii+1:size(out{i},1)
            if tolerance_matrix{ii,iii}
             % check MOB group -> possible source of errors
             if strcmp(out{i}{ii,4},out{i}{iii,4})
              % which has lower p-value      
              if out{i}{ii,6} > out{i}{iii,6}
               list_removed = [list_removed,ii];
              else
               list_removed = [list_removed,iii]; 
              end 
             end
            end
        end
    end
    % Step 3: make copy and remove elements in list
    out2{i} = out{i};
    out2{i}(list_removed,:) = [];

    %% make tmp cell
    if ~isempty(out{i})
        tmpcell = sprintf('%d,%s,%d,%s,%s,%.4f',out{i,1}{1,:});
        if size(out{i},1)>1
            for ii=2:size(out{i},1)
                tmpcell = [tmpcell,';',sprintf('%d,%s,%d,%s,%s,%.4f',...
                    out{i,1}{ii,:})];
            end
        end
    outcell{i} = tmpcell;
    end
    
    % make tmp cell out2
    if ~isempty(out2{i})
        tmpcell = sprintf('%d,%s,%d,%s,%s,%.4f',out2{i,1}{1,:});
        if size(out2{i},1)>1
            for ii=2:size(out2{i},1)
                tmpcell = [tmpcell,';',sprintf('%d,%s,%d,%s,%s,%.4f',...
                    out2{i,1}{ii,:})];
            end
        end
    outcell2{i} = tmpcell;
    end
    
    end
end

end