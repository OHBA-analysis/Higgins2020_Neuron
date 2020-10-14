function [ netmats ] = netmat_spectramt( data, S )

% data is num_nodes x num_embeddings x ntpts
% netmat is num_nodes x num_nodes x num_embeddings
% netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
   
try S.type=S.type; catch S.type='coh'; end
try S.var_normalise=S.var_normalise; catch S.var_normalise=false; end;

if isfield(S,'netmats') && isfield(S.netmats,'spectramt') && ~isempty(S.netmats.spectramt)
    netmats=S.netmats;
else

    
    Hz=S.fsample;
    
    if S.var_normalise
        data=normalise(data,2);
    end

    data=permute(data,[2 1]);

    % compute spectra
    params = struct('Fs',Hz); % Sampling rate
    params.fpass = S.fband;  % band of frequency you're interested in 
    params.tapers = [4 7]; % taper specification - leave it with default values
    %params.tapers = [5 9]; % taper specification - leave it with default values
    params.p = 0; % interval of confidence - set to 0 if you don?t wish to compute these
    %params.win = 10 * Hz; % multitaper window 
    %params.win = 5 * Hz; % multitaper window 
    %params.win = 2 * Hz; % multitaper window 
    params.win = S.reg * Hz; % multitaper window 
    
    params.to_do = [1 0]; % turn off pdc
    
    if params.win>size(data,1)
        params.win=size(data,1)-1;
    end
    
    fitmt=hmmspectramt(data,size(data,1),[],params);    
    netmats.spectramt=fitmt.state;
end;

% by default set netmat to coh, with num_embeddings in output
% corresponding to num freq bins:
netmats.netmat=permute(netmats.spectramt.(S.type),[2 3 1]);

% create netmat_full:

switch S.full_type
    case 'full' % netmat_full is (num_nodes x num_embeddings) x (num_nodes x num_embeddings)
        num_embeddings=size(netmats.netmat,3);
        num_nodes=size(netmats.netmat,1);
        netmats.netmat_full=speye(num_nodes*num_embeddings);
        for node_ind=1:num_nodes
            from = (node_ind-1)*num_embeddings+1;
            to = from+num_embeddings-1;
            for node_ind2=1:num_nodes
                from2 = (node_ind2-1)*num_embeddings+1;
                to2 = from2+num_embeddings-1;
                netmats.netmat_full(from:to,from2:to2)=sparse(diag(permute(netmats.netmat(node_ind,node_ind2,:),[3 1 2])));
            end;
        end;
    case 'mean_abs' % netmat_full is (num_nodes x num_nodes)

        netmats.netmat_full=mean(abs(netmats.netmat),3);
        
end;

netmats.type=S.type;
netmats.full_type=S.full_type;

end

