clc; clear all;

profile clear
profile('-memory','on')
profile on

%linux
addpath(genpath('C:/Users/..../SCOTSv0.2/src2'))

%Mac
addpath(genpath('/Users/..../SCOTSv0.2/src2'))

Example = 'DataMatlab';
fileID = fopen(strcat(Example,'/S_dtCon_B_labels.txt'),'r');
B_labels = fscanf(fileID,'%d');
fclose(fileID);
fileID = fopen(strcat(Example,'/S_dtCon_Bi.txt'),'r');
Bi = fscanf(fileID,'%d');
fclose(fileID);
fileID = fopen(strcat(Example,'/S_dtCon_Bj.txt'),'r');
Bj = fscanf(fileID,'%d');
fclose(fileID);
% fileID = fopen(strcat(Example,'/S_dtCon_T_row_b.txt'),'r');
% T_row_b = fscanf(fileID,'%f');
% fclose(fileID);

% T = T_row_b(1,1);
% row = T_row_b(2,1);
% b = T_row_b(3,1);

Sr = 1.3 + (1.3*1.3+20)^0.5

BN = length(B_labels);
B = sparse(Bi, Bj, ones(length(Bi),1),BN,BN);
G = digraph(B);
% plot(G)
[bins,binsizes] = conncomp(G);
Nscc = length(binsizes);    %number of scc's

% i=10
% SCCnodes = find(bins == i); %nodes in the scc
% SCC_B_labels = B_labels(SCCnodes);
% SCC_B = B(SCCnodes,SCCnodes);
% nnzero = nnz(SCC_B)  % number of nonzero elements in the sparse matrix
% SCC_G = digraph(SCC_B);
% figure(2)
% plot(SCC_G)

lev = zeros(Nscc,1);
Thn = zeros(Nscc,1);
tic;
for i=1:Nscc
    i ;
    SCCnodes = find(bins == i); %nodes in the scc
    SCC_B_labels = B_labels(SCCnodes);
    SCC_B = B(SCCnodes,SCCnodes);
    if ~(nnz(SCC_B) == 0)   % nnz(): number of nonzero elements in the sparse matrix
        SCC_dtC_Npartitions = length(unique(SCC_B_labels)); % SCC number of dtControl partitions
        SCC_dtC_partitions = unique(SCC_B_labels);
%/        comp_hgtm_dtControl takes transpose of SCC_B /%
        R = comp_hgtm_dtControl(SCC_dtC_Npartitions, SCC_B_labels, SCC_B', SCC_dtC_partitions)
        lev(i,1) = eigs(R,1);   %eigenvalue with largest magnitude
        if isnan(lev(i,1))
            lev(i,1) = max(abs(eigs(R)));
            if isnan(lev(i,1))
                lev(i,1) = max(abs(eig(full(R))));
            end
        end
        hn(i,1) = log2(lev(i,1));
    end
end
toc
% lev
% hn
dtcNpartitions = length(unique(B_labels))
Thn = log2(max(lev))

profile viewer

% 
% 
% dtC_partitions = unique(B_labels);
% dtC_Npartitions = length(unique(B_labels));
% R = comp_hgtm_dtControl(dtC_Npartitions, B_labels, B, dtC_partitions)


