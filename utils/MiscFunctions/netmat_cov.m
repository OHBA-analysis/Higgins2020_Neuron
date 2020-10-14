function [ netmats ] = netmat_cov( data, S )

% data is num_nodes x num_embeddings x ntpts
% netmats.netmat is num_nodes x num_nodes x num_embeddings
% netmats.netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
%       where the (num_nodes x num_embeddings) dimension is in the order:
%       [node1,embedding1; node1,embedding2 ... node1,embeddingN;
%       [node2,embedding1; node2,embedding2 ... node2,embeddingN;
%       etc...

try S.type=S.type; catch S.type='full'; end; % 'full' or 'partial'
try S.var_normalise=S.var_normalise; catch S.var_normalise=false; end;

% remove any power/variance info
% but note that this does not remove autocorr info:

if S.var_normalise
    data=normalise(data,3);
end;

netmats.netmat = zeros(size(data,1),size(data,1),size(data,2));

% compute netmat separately for each embedding
for kk=1:size(data,2),
    dat=permute(data(:,kk,:),[1 3 2]);
    netmats.netmat(:,:,kk)=dat*dat'/size(dat,2);
end;

% compute netmat combined over all embeddings
if size(data,2)>1
    dat=reshape(permute(data,[2 1 3]),[size(data,1)*size(data,2),size(data,3)]);
    netmats.netmat_full=dat*dat'/size(dat,2);
else
    netmats.netmat_full=netmats.netmat(:,:,1);    
end;

if strcmp(S.type,'partial')
    netmats.netmat_full=pinv(netmats.netmat_full);
end;

end

