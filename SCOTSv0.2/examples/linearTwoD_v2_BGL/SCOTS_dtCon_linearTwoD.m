clc; clearvars;

profile clear
profile('-memory','on')
profile on

addpath(genpath('../../mfiles/mexfiles'))

Example = 'DataMatlab';
fileID = fopen(strcat(Example,'/S_dtCon_B_NPosts.txt'),'r');
B_NPosts = fscanf(fileID,'%d');
fclose(fileID);
fileID = fopen(strcat(Example,'/Xi.txt'),'r');
Xi = fscanf(fileID,'%d');
fclose(fileID);
fileID = fopen(strcat(Example,'/Xj.txt'),'r');
Xj = fscanf(fileID,'%d');
fclose(fileID);

BN = length(B_NPosts);

X = unique(Xi);
Bi = zeros(length(Xi),1);
Bj = zeros(length(Xj),1);
for k=1:length(X)
    Bi(Xi==X(k)) = k;
    Bj(Xj==X(k)) = k;
end

B = sparse(Bi, Bj, ones(length(Bi),1),BN,BN);
G = digraph(B);
% plot(G)
[bins,binsizes] = conncomp(G);
Nscc = length(binsizes);    %number of scc's
max_scc_i_NPosts = zeros(Nscc,1);
for k=1:Nscc
    SCCnodes = find(bins == k); %indices of the nodes in the scc
    SCC_B_NPosts = B_NPosts(SCCnodes);
    max_scc_i_NPosts(k) = max(SCC_B_NPosts);
    %SCC_B = B(SCCnodes,SCCnodes);
end
wmax = log2(max(max_scc_i_NPosts))
    


%profile viewer

%%
% % i=10
% % SCCnodes = find(bins == i); %nodes in the scc
% % SCC_B_labels = B_labels(SCCnodes);
% % SCC_B = B(SCCnodes,SCCnodes);
% % nnzero = nnz(SCC_B)  % number of nonzero elements in the sparse matrix
% % SCC_G = digraph(SCC_B);
% % figure(2)
% % plot(SCC_G)
% 
% lev = zeros(Nscc,1);
% Thn = zeros(Nscc,1);
% tic;
% for i=1:Nscc
%     i ;
%     SCCnodes = find(bins == i); %nodes in the scc
%     SCC_B_labels = B_labels(SCCnodes);
%     SCC_B = B(SCCnodes,SCCnodes);
%     if ~(nnz(SCC_B) == 0)   % nnz(): number of nonzero elements in the sparse matrix
%         SCC_dtC_Npartitions = length(unique(SCC_B_labels)); % SCC number of dtControl partitions
%         SCC_dtC_partitions = unique(SCC_B_labels);
% %/        comp_hgtm_dtControl takes transpose of SCC_B /%
%         R = comp_hgtm_dtControl(SCC_dtC_Npartitions, SCC_B_labels, SCC_B', SCC_dtC_partitions)
%         lev(i,1) = eigs(R,1);   %eigenvalue with largest magnitude
%         if isnan(lev(i,1))
%             lev(i,1) = max(abs(eigs(R)));
%             if isnan(lev(i,1))
%                 lev(i,1) = max(abs(eig(full(R))));
%             end
%         end
%         hn(i,1) = log2(lev(i,1));
%     end
% end
% toc
% % lev
% % hn
% dtcNpartitions = length(unique(B_labels))
% Thn = log2(max(lev))
%%

%profile viewer


