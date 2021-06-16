clc; clearvars;

profile clear
profile('-memory','on')
profile on

%win10
addpath(genpath('C:/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/src2'))

%Mac
addpath(genpath('/Users/mst/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/src2'))

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


