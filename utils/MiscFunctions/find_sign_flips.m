function [ sign_flipped_spm_fnames, sign_flip_results ] = find_sign_flips( Sin )

% [ sign_flipped_spm_fnames, sign_flip_results ] = find_sign_flips( S )
%
% S.Ds=Ds
% S.num_iters=1000;
% S.roinets_protocol='symmetric'
% S.prefix='sf_'

Ds=Sin.Ds;

try num_uber_updates=Sin.num_uber_updates; catch num_uber_updates=2*length(Ds); end;
try method=Sin.method; catch method='matrix_distance_metric'; end;
try roinets_protocol=Sin.roinets_protocol; catch roinets_protocol='symmetric'; end;
try innovations_mar_order=Sin.innovations_mar_order; catch innovations_mar_order=14; end;
try num_iters_within_updates=Sin.num_iters; catch num_iters=100; end;
try num_embeddings=Sin.num_embeddings; catch num_embeddings=10; end;
try subj_template=Sin.subj_template; catch subj_template=1; end;

%% compute embedded cov matrices
S=[];
S.concat = [];
S.concat.protocol=roinets_protocol;
S.concat.innovations_mar_order = innovations_mar_order;
S.concat.embed.do=1;
S.concat.embed.num_embeddings=num_embeddings;
S.concat.embed.rectify=false;
S.concat.whiten=1;
S.concat.normalisation='voxelwise';
S.concat.pcadim=-1;
S.netmat_method=@netmat_cov;

[ state_netmats_cov ] = hmm_full_global_cov( Ds, S );

state_netmats_cov_in=state_netmats_cov;

num_embeddings=size(state_netmats_cov{1}.global.netmat,3);
num_nodes=size(state_netmats_cov{1}.global.netmat,1);
num_lags=num_embeddings;

clear netmats;
for subj=1:length(state_netmats_cov),
    netmats(:,:,subj)=state_netmats_cov{subj}.global.netmat_full;
end;

% netmats
num_subj=size(netmats,3);

flips=ones(num_nodes,num_subj);
   
% initialise
max_num_nodes_to_change=20;
naccept=0;
nreject=0;
sigma=zeros(num_subj,1);
energies=[];
old_energies=[];
new_energies=[];
rejections=[];

flips=ones(num_nodes,num_subj);
best_flips=ones(num_nodes,num_subj);
best_energy=ones(num_subj,1)*1e32;
sub_jump_count=zeros(num_subj,1);

% initialise global cov to first subject
group_netmat=netmats(:,:,subj_template);

uber_mat=netmats;
stick=false;

for uit=1:num_uber_updates
      
    uit/num_uber_updates
    
    pinv_group_netmat=pinv(group_netmat);
    inds=[];
    for subj=1:num_subj,
        
        if subj~=subj_template

            [old_energy, inds]=sign_flip_post(netmats(:,:,subj),flips(:,subj),group_netmat,sigma(subj),num_embeddings,pinv_group_netmat, method,num_lags,inds);

            for sit=1:num_iters_within_updates

                % propose new flips
                flips_old=flips(:,subj);

                if sit==1 && uit==1
                    num_nodes_to_change=num_nodes; % for first try flip everything
                else
                    num_nodes_to_change=randperm(max_num_nodes_to_change,1);            
                end;
                nodeinds=randperm(num_nodes,num_nodes_to_change);            
                flips(nodeinds,subj)=-flips(nodeinds,subj);

                new_energy=sign_flip_post(netmats(:,:,subj),flips(:,subj),group_netmat,sigma(subj),num_embeddings,pinv_group_netmat, method,num_lags,inds);

                %a=exp((old_energy-new_energy)/temperature);            

                reject=true;
                %if a > rand(1)	% Accept the new state.
                if new_energy<old_energy
                %a
                    %old_energy-new_energy
                    %old_energy
                    %new_energy
                    reject=false;
                end;


                sub_jump_count(subj)=sub_jump_count(subj)+1;
                old_energies(sub_jump_count(subj),subj)=old_energy;
                new_energies(sub_jump_count(subj),subj)=new_energy;
                rejections(sub_jump_count(subj),subj)=reject;

                if reject
    %                disp('reject')

                    nreject=nreject+1;

                    flips(:,subj)=flips_old;                

                else
                    old_energy=new_energy;
                    naccept=naccept+1;

                    if new_energy<best_energy(subj)
                        best_flips(:,subj)=flips(:,subj);
                        best_energy(subj)=new_energy;
                    end;
                end;

                energies(sub_jump_count(subj),subj)=old_energy;

            end; % for sit=1:num_iters_within_updates  
        end;
    end; % for subj=1:num_subj,
    
%     % update group netmat
%     for subj=1:num_subj,
%         flips_subj=sparse(diag(kron(best_flips(:,subj)',ones(1,num_embeddings))));
%         uber_mat(:,:,subj)=flips_subj*netmats(:,:,subj)*flips_subj;
%     end;
%      
%     if uit<num_subj,
%         subj_template=subj_template+1;
%     else
%         if ~stick
%             % find best template and stick with it
%             ens=sum(energies(num_iters_within_updates:num_iters_within_updates:end,:),2);
% 
%             [val subj_template]=min(ens);
%             
%             % reset the best energies:
%             best_energy=ones(num_subj,1)*1e32;
% 
%             stick=true;
%         end;
%             
%     end;
%     
%     group_netmat=netmats(:,:,subj_template);

    disp(['Using subj ' num2str(subj_template) ' for template']);
    
    % update proposal based on rejection rate
    rej_rate=nreject/(nreject+naccept);
    %disp(['rej_rate=' num2str(rej_rate)]);
    
    %temperature=temperature+round(temperature*(rej_rate-0.6));
    
    naccept=0;
    nreject=0;
    
    %disp(['temperature=' num2str(temperature)]);

end;

%% apply sign flipping to prepared parcellated files

nsubs=length(Ds);

for subnum = 1:nsubs       
    Dp = spm_eeg_load(Ds{subnum});
    
    Dsnew{subnum}=Dp;
    
    % create new object
    S=[];
    S.D=Dsnew{subnum};
    S.outfile = prefix(fullfile(Ds{subnum}),Sin.prefix);
    Dsnew{subnum}=spm_eeg_copy(S);
    
    sign_flipped_spm_fnames{subnum}=Dsnew{subnum}.fullfile;
    
    for pp=1:size(Dp,1)
        for tri=1:size(Dsnew{subnum},3)
            Dsnew{subnum}(pp,:,tri)=Dp(pp,:,tri)*best_flips(pp,subnum);            
        end;       
    end;
    
    Dsnew{subnum}.save;
end

%% compute cov after flipping
S=[];

S.concat = [];
S.concat.protocol=roinets_protocol;
S.concat.innovations_mar_order = innovations_mar_order;
S.concat.embed.do=1;
S.concat.embed.num_embeddings=num_embeddings;
S.concat.embed.rectify=false;
S.concat.whiten=1;
S.concat.normalisation='voxelwise';
S.concat.pcadim=-1;
S.netmat_method=@netmat_cov;

[ state_netmats_cov_out ] = hmm_full_global_cov( Dsnew, S );

if 0,
    %
    
    for subnum = 1:nsubs     
        flipstmp=best_flips(:,subnum);
        flips=sparse(diag(kron(flipstmp',ones(1,num_embeddings))));
        state_netmats_cov_flipped2{subnum}.global.netmat_full=flips*state_netmats_cov{subnum}.global.netmat_full*flips;
    end;
    
    %
    
    figure;imagesc(triu((state_netmats_cov_flipped{2}.global.netmat_full),15))
    figure;imagesc(triu((state_netmats_cov_flipped2{2}.global.netmat_full),15))

    
    matrix_distance_metric(state_netmats_cov_flipped{2}.global.netmat_full, state_netmats_cov_flipped{1}.global.netmat_full,15,'none')

end;

sign_flip_results.state_netmats_cov_in=state_netmats_cov_in;
sign_flip_results.state_netmats_cov_out=state_netmats_cov_out;
sign_flip_results.best_flips=best_flips; 
sign_flip_results.energies=energies;
sign_flip_results.sign_flipped_spm_fnames=sign_flipped_spm_fnames;
