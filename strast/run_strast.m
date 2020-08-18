%% import data

cd '../data'
fn_fasta = dir('*.fasta');
plasmids = cell(length(fn_fasta),2);
targetSet = cell(length(fn_fasta),1);

for i=1:length(fn_fasta)
   plasmids{i,1} = fastaread(fn_fasta(i).name);
   plasmids{i,1}.Sequence = upper(plasmids{i,1}.Sequence);
   
   targetSet{i} = plasmids{i,1}.Sequence;
   disp(length(targetSet{i}))
end

%% run strast

cd '../strast'

% set params
vals = [4,6];
vals2 = [4,6];
lnic = [140];
tresh= -[12,20,30];

% load parameters
load ../data/dataset_nic_112full_1_11_17.mat nic112pval
querySet = nic112pval;
load ../data/idx_cdist_MAT_26_10_17.mat CdistMAT idxMAT
load ../data/p_poly_C6_6_11_17.mat p_poly_C6

name=sprintf('q112_t1_s7c128');
[out,~] = strast_mainV6_2(targetSet,querySet,idxMAT,CdistMAT,vals,p_poly_C6,vals2,tresh,lnic);
eval([name,'=out;']);
clear out

%% parse results

load ../data/dataset_nic_112full_1_11_17.mat NIC112
querytable = NIC112;
[~,~,~,outcell2] = strastV1OutputTidyV4_1(q112_t1_s7c128,tresh,querytable);

% store results
for i=1:length(outcell2)
    outcell2{i} = [fn_fasta(i).name,';',outcell2{i}];
end
cell2csv('../data/predictions.csv',outcell2)
