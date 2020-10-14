function hfig=plot_surf_summary_neuron(p,data,surface_inflation,single_plot,interptype,view_angle,thresh,clims)
    % Render volume data onto a surface
    %
    % *REQUIRES WORKBENCH*
    %
    % INPUTS
    % - p - parcellation object
    % - data - data that can be saved to a nii file via p.savenii(data)
    % - surface_inflation - integer level of inflation for display surface (default=0, no inflation)
    % - single_plot - put the two hemispheres together (default=false, like Workbench)
    % - interptype - passed to workbench (default='trilinear')
    % - view_angle - passed to view to set rotation

    if nargin < 8
        clims = [];
    end
    
    if nargin < 7 || isempty(thresh)
        thresh = 0;
    end
    
    if nargin < 6 || isempty(view_angle)
        view_angle = [-90 0];
    end

    if nargin < 5 || isempty(interptype)
        interptype = 'trilinear';
    end

    if nargin < 4 || isempty(single_plot)
        single_plot = false;
    end

    if nargin < 3 || isempty(surface_inflation)
        surface_inflation = 0;
    end

    if single_plot == true && surface_inflation ~= 0
        fprintf(2,'Warning, single plot with inflated surface does not currently plot correctly');
    end

    niifile = p.savenii(data);
    output_right    = [niifile '_right.func.gii'];
    output_left     = [niifile '_left.func.gii'];
    cl = onCleanup(@()  cellfun(@delete,{niifile,output_left,output_right})); % Enable deleting temp files even if debugging

    surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii');
    surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii');

    switch surface_inflation
        case 0
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii');
        case 1
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.inflated.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.inflated.32k_fs_LR.surf.gii');
        case 2
            display_surf_right = fullfile(osldir,'std_masks','ParcellationPilot.R.very_inflated.32k_fs_LR.surf.gii');
            display_surf_left = fullfile(osldir,'std_masks','ParcellationPilot.L.very_inflated.32k_fs_LR.surf.gii');
    end

    % Map volume to surface
%     runcmd('/Applications/workbench/bin_macosx64/wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_right,output_right,interptype)
%     runcmd('/Applications/workbench/bin_macosx64/wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_left,output_left,interptype)
     runcmd('/Users/chiggins/.local/src/workbench/bin_macosx64/wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_right,output_right,interptype)
     runcmd('/Users/chiggins/.local/src/workbench/bin_macosx64/wb_command -volume-to-surface-mapping %s %s %s -%s',niifile,surf_left,output_left,interptype)

    sl = gifti(display_surf_left);
    vl = gifti(output_left);
    sr = gifti(display_surf_right);
    vr = gifti(output_right);
    hfig = figure('Position',[100 100 1024 256]);
    set(hfig,'Color','White');
    if isempty(clims)
        clims = [min([min(vl.cdata) min(vr.cdata)]) max([max(vl.cdata) max(vr.cdata)])];
    end
    vl.cdata(vl.cdata==0)=NaN;
    vr.cdata(vr.cdata==0)=NaN;
    % make colormap
    if min([min(vl.cdata) min(vr.cdata)]) < 0 || clims(1) < 0
        cm = cat(2,linspace(.5, 1 ,63)',linspace(0, 1 ,63)',linspace(0, 0 ,63)');
        cm2 = cat(2,linspace(0, 0 ,63)',linspace(1, 0 ,63)',linspace(1, .5 ,63)');
        %cm = cat(1,cm2,[.6 .6 .6],cm);
        cm = cat(1,flipud(cm2),flipud(cm));
    else
        cm = cat(2,linspace(.5, 1 ,63)',linspace(0, 1 ,63)',linspace(0, 0 ,63)');
        c=flipud(cm);
        cm = cat(1,[.6 .6 .6],cm);
    end
    
    for ii = 1:5
        if ii==3
            ax(ii) = axes('Position',[.05+((ii-1)*.165) .1 .175 .8]);
        elseif ii==2
                
            ax(ii) = axes('Position',[.08+((ii-1)*.165) .3 .175 .5]);
        elseif ii==5
            ax(ii) = axes('Position',[.025+((ii-1)*.165) .3 .175 .5]);
        else
            if ii==1 || surface_inflation==0
                ax(ii) = axes('Position',[.05+((ii-1)*.165) .3 .175 .5]);
            else
                ax(ii) = axes('Position',[.05+((ii-1)*.165)-0.05 .3 .175 .5]);
            end
        end
    end
    
    % lateral
    axes(ax(1));
    s(1) = patch('Faces',sl.faces,'vertices',sl.vertices,'CData',[]);
    hold on
    sg(1) = patch('Faces',sl.faces,'vertices',sl.vertices);
    
    axes(ax(2));
    s(2) = patch('Faces',sr.faces,'vertices',sr.vertices,'CData',[]);
    hold on
    sg(2) = patch('Faces',sr.faces,'vertices',sr.vertices);    
    
    set(s(1),'FaceVertexCData',vl.cdata)
    set(s(2),'FaceVertexCData',vr.cdata)

    set(sg(1),'FaceVertexCData',0.4*ones(size(vl.cdata,1),3));
    set(sg(2),'FaceVertexCData',0.4*ones(size(vr.cdata,1),3));
    set(sg(1),'FaceVertexAlphaData',+~isfinite(vl.cdata),'FaceAlpha','interp','AlphaDataMapping','none');
    set(sg(2),'FaceVertexAlphaData',+~isfinite(vr.cdata),'FaceAlpha','interp','AlphaDataMapping','none');    
    
    view(ax(1),[270 0])
    view(ax(2),[-270 0])
    clear s sg

    % top view
    if surface_inflation==0
            
    axes(ax(3));
    s(1) = patch('Faces',sl.faces,'vertices',sl.vertices,'CData',[]);
    hold on
    sg(1) = patch('Faces',sl.faces,'vertices',sl.vertices);
    s(2) = patch('Faces',sr.faces,'vertices',sr.vertices,'CData',[]);
    sg(2) = patch('Faces',sr.faces,'vertices',sr.vertices);
    
    set(s(1),'FaceVertexCData',vl.cdata)
    set(s(2),'FaceVertexCData',vr.cdata)

    set(sg(1),'FaceVertexCData',0.4*ones(size(vl.cdata,1),3));
    set(sg(2),'FaceVertexCData',0.4*ones(size(vr.cdata,1),3));
    set(sg(1),'FaceVertexAlphaData',+~isfinite(vl.cdata),'FaceAlpha','interp','AlphaDataMapping','none');
    set(sg(2),'FaceVertexAlphaData',+~isfinite(vr.cdata),'FaceAlpha','interp','AlphaDataMapping','none');    
    end
    clear s sg

    % medial
    axes(ax(4));
    s(1) = patch('Faces',sl.faces,'vertices',sl.vertices,'CData',[]);
    hold on
    sg(1) = patch('Faces',sl.faces,'vertices',sl.vertices);
    
    axes(ax(5));
    s(2) = patch('Faces',sr.faces,'vertices',sr.vertices,'CData',[]);
    hold on
    sg(2) = patch('Faces',sr.faces,'vertices',sr.vertices);    
    
    set(s(1),'FaceVertexCData',vl.cdata)
    set(s(2),'FaceVertexCData',vr.cdata)

    set(sg(1),'FaceVertexCData',0.4*ones(size(vl.cdata,1),3));
    set(sg(2),'FaceVertexCData',0.4*ones(size(vr.cdata,1),3));
    set(sg(1),'FaceVertexAlphaData',+~isfinite(vl.cdata),'FaceAlpha','interp','AlphaDataMapping','none');
    set(sg(2),'FaceVertexAlphaData',+~isfinite(vr.cdata),'FaceAlpha','interp','AlphaDataMapping','none');    
    
    view(ax(4),[90 0])
    view(ax(5),[-90 0])
    clear s sg
    
    
    arrayfun(@(x) shading(x,'interp'),ax);
    arrayfun(@(x) axis(x,'off'), ax);

    
    for ii = 1:5
        axes(ax(ii))
        colormap(cm)
        %zoom(1.1)
        if ~isempty(clims)
            caxis(clims);
        else
            mn = min([min(vl.cdata) min(vr.cdata)]);
            mx = max([max(vl.cdata) max(vr.cdata)]);
            if mn < 0
                m = max([ abs(mn) abs(mx) ]);
                caxis([-m m]);
            else
                caxis([0 mx])
            end
        end
        material( [0.3, 0.8, 0.2] );
        camlight('left')
        camlight('right')
        camlight('headlight')
        if ii == 5
            c = colorbar('EastOutside');
            c.Position(1) = .9;
            c.Position(2) = .25;
            c.Position(3) = .02;
            c.FontSize = 18;
        end
    end
    
    
end