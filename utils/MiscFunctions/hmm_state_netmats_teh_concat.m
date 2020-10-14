function [ state_netmats ] = hmm_state_netmats_concat( hmm, S )

% [ state_netmats ] = hmm_state_netmats( hmm, S )

try 
    S.global_only=S.global_only; 
catch
    S.global_only=false;
end

logtrans=0;
state_netmats=[];

normalisation = S.normalisation;
embed=S.embed;
assignment=S.assignment;
subj_inds=hmm.subj_inds;

D_fnames=S.parcellated_filenames;

nsubs=length(D_fnames);

datap=[];

for subnum = 1:nsubs

    disp(['Concatenating for subj num ' num2str(subnum)]);
    try
        D = spm_eeg_load(prefix(D_fnames{subnum},'p'));
    catch
        D = spm_eeg_load((D_fnames{subnum}));
    end
    
    if subnum==1
        state_netmats{1}.parcellation=D.parcellation;
    end

    embed.tres=1/D.fsample;

    %data = prepare_data(D,normalisation,logtrans,embed,roinets_protocol);
    data = osl_teh_prepare_data(D,normalisation,logtrans,[],embed);
    
    %data=permute(data,[1 3 2]);
    
    datap = [datap, data];
  
end

if ~S.global_only && hmm.K>1
        
    hmm_sub = hmm; 
    hmm_sub.statepath = hmm.statepath(subj_inds<=nsubs); 
    hmm_sub.gamma = hmm.gamma(subj_inds<=nsubs,:);

    for k = 1:hmm_sub.K

        disp(['Computing for state ' num2str(k)]);

        switch assignment
        case 'hard' 

            inds = logical(hmm_sub.statepath == k);
            
            if size(datap,2) ~= length(inds)
                error('not implemented');
            end
            
            datap_in=datap(:,inds);
                          
         case 'soft'
            x = hmm_sub.gamma(:,k);

            if size(datap,2) ~= length(x)
                error('not implemented');
            end

            x2=permute(repmat(x,[1, size(datap,1)]),[2 1]);
            datap_in=datap.*x2;
        end

        [state_netmats{1}.state{k}]=feval(S.netmat_method,datap_in,S.netmat_method_options);

    end

end
    
[state_netmats{1}.global]=feval(S.netmat_method,datap,S.netmat_method_options);

%sum(state_netmats{subnum}.global.spectramt.psd(:,9,9),1)
%sum(state_netmats{subnum}.state{k}.spectramt.psd(:,9,9),1)
%std(squeeze(datap_in(9,1,:))),std(squeeze(datap(9,1,:)))
%sum(state_netmats{subnum}.global.netmat(9,9,:))
%sum(state_netmats{subnum}.state{k}.netmat(9,9,:))


state_netmats{1}.netmat_method=S.netmat_method;
                        
