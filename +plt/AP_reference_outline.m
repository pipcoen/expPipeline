function h = AP_reference_outline(type,color,reference_im)
% h = AP_reference_outline(type,color,reference_im)
%
% If a reference image is provided, bregma is defined there
% type - grid, ccf, ccf_wf, ccf_aligned, retinotopy (master only at the moment)

%% Define bregma
um2pixel = 20.6;
    case 'ccf'
        % Plot allen CCF outlines       
        load('allenCorticalBoundaries.mat');
        bregma = allenCCFbregma;
        for curr_area_idx =1:length(cortical_area_boundaries)
            h = cellfun(@(outline) plot((outline(:,2)-bregma(3))*10, ...
                (bregma(1)-outline(:,1))*10,'color',color),cortical_area_boundaries{curr_area_idx},'uni',false);
        end       
end


% ONLY NEEDED ONCE: create top-down cortical boundaries from CCF
% % Get first brain pixel from top-down, get annotation at that point
% [~,top_down_depth] = max(av>1, [], 2);
% top_down_depth = squeeze(top_down_depth);
%
% [xx,yy] = meshgrid(1:size(top_down_depth,2), 1:size(top_down_depth,1));
% top_down_annotation = reshape(av(sub2ind(size(av),yy(:),top_down_depth(:),xx(:))), size(av,1), size(av,3));
%
% % Get all labelled areas
% used_areas = unique(top_down_annotation(:));
%
% % Restrict to only cortical areas
% structure_id_path = cellfun(@(x) textscan(x(2:end),'%d', 'delimiter',{'/'}),st.structure_id_path);
%
% ctx_path = [997,8,567,688,695,315];
% ctx_idx = find(cellfun(@(id) length(id) > length(ctx_path) & ...
%     all(id(min(length(id),length(ctx_path))) == ctx_path(min(length(id),length(ctx_path)))),structure_id_path));
%
% plot_areas = intersect(used_areas,ctx_idx);
%
% bregma = allenCCFbregma;
%
% % Get outlines of all areas
% cortical_area_boundaries = cell(size(plot_areas));
% for curr_area_idx = 1:length(plot_areas)
%     cortical_area_boundaries{curr_area_idx} = bwboundaries(top_down_annotation == plot_areas(curr_area_idx));
% end













