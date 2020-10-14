function [ netmats ] = hmm_full_global_cov( Ds, S )

logtrans=0;

normalisation  = S.concat.normalisation;
roinets_protocol=S.concat.protocol;
embed=S.concat.embed;
try innovations_mar_order=S.innovations_mar_order; catch, innovations_mar_order=14; end;

nsubs=length(Ds);

clear netmats;

for subnum = 1:nsubs

    disp(['Computing for subj num ' num2str(subnum)]);
       
    Dp = spm_eeg_load(Ds{subnum});
        
    embed.tres=1/Dp.fsample;
            
    % returns data as num_nodes x num_embeddings x ntpts
    [datap, num_embeddings] = prepare_data(Dp,normalisation,logtrans,embed,roinets_protocol,innovations_mar_order);        
       
    netmats{subnum}.global=netmat_cov( datap, S );
    
    netmats{subnum}.num_embeddings=num_embeddings;
end                        
