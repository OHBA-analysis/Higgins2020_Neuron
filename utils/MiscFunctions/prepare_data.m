function [data, num_embeddings] = prepare_data(D,normalisation,logtrans,embed,roinets_protocol,innovations_mar_order)

if nargin < 6
    innovations_mar_order=14;
end

% Returns data as num_nodes x num_embeddings x ntpts

% Reshape trialwise data
data = D(:,:,:);
data = reshape(data,[D.nchannels,D.nsamples*D.ntrials]);

% Select only good data
good_samples = ~all(badsamples(D,':',':',':'));
good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);
data = data(:,good_samples);

% Log transform
if logtrans
    data = log10(data);
end

% Do orthogonalisation
if ~strcmp(roinets_protocol,'none'),
    if ~strcmp(roinets_protocol, 'innovations_mar')
        data = ROInets.remove_source_leakage(data,roinets_protocol);
    else            
        data = leakcorr(data',size(data,2),innovations_mar_order)';
    end
end

% Do embedding
if exist('embed') && embed.do,    
    
    disp('Time embedding data');
    
    if isfield(embed,'centre_freq')
        span=1/embed.centre_freq; %secs    
        num_embeddings=round(span/embed.tres); 
    elseif isfield(embed,'num_embeddings')
        num_embeddings=embed.num_embeddings; 
    end

    lags=round(linspace(-num_embeddings/2,num_embeddings/2,num_embeddings));
    %lags=round(-num_embeddings/2:num_embeddings/2);
    %lags=round(0:num_embeddings-1);

    disp(['lags=' mat2str(lags)]);
    dataout=randn(size(data,1),num_embeddings,size(data,2))*std(squash(data(:,1:500))); % num_nodes x num_embeddings x ntpts
    for k=1:size(data,1)
        [tmp,valid]=embedx(data(k,:)',lags);
        dataout(k,:,valid)=tmp';
    end;
      
    data=dataout;
    clear dataout;
else
    data=permute(data,[1 3 2]);
    num_embeddings=1;
end;

% Normalisation
switch normalisation
    case 'global'
        data = demean(data,3)./std(data(:));
    case 'voxelwise'
        data = normalise(data,3);
    case 'none'
        data = demean(data,3);
end

if embed.rectify
    data=abs(data);
    
    switch normalisation
    case 'global'
        data = demean(data,3)./std(data(:));
    case 'voxelwise'
        data = normalise(data,3);
    case 'none'
        data = demean(data,3);
    end;
end;

end

