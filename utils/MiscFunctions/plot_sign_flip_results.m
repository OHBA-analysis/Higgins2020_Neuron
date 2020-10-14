function plot_sign_flip_results(state_netmats_cov_glob,state_netmats_cov_glob_flip, template_subj, freq, sign_flip_results, print_fname)

if nargin < 4
    freq=[];
end

if nargin < 5
    sign_flip_results=[];
end

if nargin < 5
    print_fname = [];
end

%% compare matrices between subjects pre-flipping

diag_offset = state_netmats_cov_glob{1}.num_embeddings;

do_post = ~isempty(state_netmats_cov_glob_flip);

state_netmats=state_netmats_cov_glob;

modes={'none','abs'};

metric_global=zeros(length(state_netmats),length(state_netmats),length(modes));

for mm=1:length(modes),
    for subj=1:length(state_netmats)
       for subj2=1:length(state_netmats)
            if subj2~=subj
                metric_global(subj,subj2,mm)=matrix_distance_metric(state_netmats{subj}.global.netmat_full, state_netmats{subj2}.global.netmat_full,diag_offset,modes{mm},[]);
            end;
       end; 
    end;
end;

figure; 
set(gcf,'Position',[440   378   1300   420]);
    
subplot_grid(2,4,1,1);
imagesc(metric_global(:,:,1),[0 0.3]);colorbar;
if ~isempty(freq)
    title(['Preflip raw: ' num2str(freq(1)) ' to ' num2str(freq(2)) 'hz']);
end

subplot_grid(2,4,1,2);
imagesc(metric_global(:,:,2),[0 0.6]);colorbar;
if ~isempty(freq)
    title(['Preflip env: ' num2str(freq(1)) ' to ' num2str(freq(2)) 'hz']);
end

% compare matrices between subjects post-flipping
if do_post

    state_netmats=state_netmats_cov_glob_flip;

    metric_global_post=zeros(length(state_netmats),length(state_netmats),length(modes));

    diag_offset=15;

    for mm=1:length(modes),
        for subj=1:length(state_netmats)
           for subj2=1:length(state_netmats)
                if subj2~=subj
                    metric_global_post(subj,subj2,mm)=matrix_distance_metric(state_netmats{subj}.global.netmat_full, state_netmats{subj2}.global.netmat_full,diag_offset,modes{mm},[]);
                end;
           end; 
        end;
    end;

    subplot_grid(2,4,2,1);
    imagesc(metric_global_post(:,:,1),[0 0.3]);colorbar;    
    if ~isempty(freq)
        title(['Postflip raw: ' num2str(freq(1)) ' to ' num2str(freq(2)) 'hz']);
    end
    
    subplot_grid(2,4,2,2);
    imagesc(metric_global_post(:,:,2),[0 0.6]);colorbar;
    if ~isempty(freq)
        title(['Postflip env: ' num2str(freq(1)) ' to ' num2str(freq(2)) 'hz']);
    end
end;

% plot covs for template subject before flipping
subplot_grid(2,4,1,3);
x=state_netmats_cov_glob{template_subj}.global.netmat_full;
x = triu(x,diag_offset)+tril(x,-diag_offset);
imagesc(x);colorbar;
title(['Template subj cov mat, preflip']);

% plot covs for template subject after flipping
subplot_grid(2,4,2,3);
x=state_netmats_cov_glob_flip{template_subj}.global.netmat_full;
x = triu(x,diag_offset)+tril(x,-diag_offset);
imagesc(x);colorbar;
title(['Template subj cov mat, postflip']);

% plot differences in covs to template subject
subplot_grid(2,4,1,4);
plot(metric_global(template_subj,1:end,1),'r','LineWidth',2); ho; 
if do_post
    plot(metric_global_post(template_subj,1:end,1),'g','LineWidth',2);
end;

plot(metric_global(template_subj,1:end,2),'b','LineWidth',2);
if do_post
    legend('pre-flip','post-flip','env');
else
    legend('pre-flip','env');    
end
plot4paper('subj #','Distance Metric');
if ~isempty(freq)
    title([num2str(freq(1)) ' to ' num2str(freq(2)) 'hz']);
end

% plot convergence
if ~isempty(sign_flip_results)
    subplot_grid(4,4,3,4);
    plot(sign_flip_results.energies,'LineWidth',2);
    plot4paper('iterations','metric');

    subplot_grid(4,4,4,4);
    plot(sign_flip_results.energies_group,'LineWidth',2);
    plot4paper('iterations','metric');
end

if ~isempty(print_fname)
    print(print_fname,'-dpng');
end

