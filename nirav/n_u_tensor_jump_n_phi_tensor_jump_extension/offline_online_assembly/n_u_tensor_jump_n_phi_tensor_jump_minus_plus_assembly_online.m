function [ res ] = ...
    n_u_tensor_jump_n_phi_tensor_jump_minus_plus_assembly_online...
    ( offline_assembly_minus_plus, para_mapping, params, grid, el_subd)

gids = ldg_global_dof_index(params, grid);

minus_plus_11_assembly = offline_assembly_minus_plus.minus_plus_11;
minus_plus_12_assembly = offline_assembly_minus_plus.minus_plus_12;
minus_plus_21_assembly = offline_assembly_minus_plus.minus_plus_21;
minus_plus_22_assembly = offline_assembly_minus_plus.minus_plus_22;

% for i = 1:1:length(el_subd)
%     para_mapping{i} = poldecomp((para_mapping{i})');
% end

for i = 1:1:length(el_subd)
    para_mapping{i} = eye(2);
end

for i = 1:1:length(el_subd)
    for j = 1:1:length(el_subd)
        minus_plus_11_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) = ...
            (para_mapping{j}(1,:) * (para_mapping{i}(1,:))') * ...
            minus_plus_11_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) / ...
            (det(para_mapping{i}) * det(para_mapping{j}));
        minus_plus_12_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) = ...
            (para_mapping{j}(1,:) * (para_mapping{i}(2,:))') * ...
            minus_plus_12_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) / ...
            (det(para_mapping{i}) * det(para_mapping{j}));
        minus_plus_21_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) = ...
            (para_mapping{j}(2,:) * (para_mapping{i}(1,:))') * ...
            minus_plus_21_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) / ...
            (det(para_mapping{i}) * det(para_mapping{j}));
        minus_plus_22_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) = ...
            (para_mapping{j}(2,:) * (para_mapping{i}(2,:))') * ...
            minus_plus_22_assembly(gids(el_subd{i},:),gids(el_subd{j},:)) / ...
            (det(para_mapping{i}) * det(para_mapping{j}));
    end
end

res = minus_plus_11_assembly + minus_plus_12_assembly + ...
    minus_plus_21_assembly + minus_plus_22_assembly;

end