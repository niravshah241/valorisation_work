function [ res ] = ...
    n_u_tensor_jump_n_phi_tensor_jump_plus_plus_assembly_online...
    ( offline_assembly_plus_plus, para_mapping, params, grid, el_subd)

gids = ldg_global_dof_index(params, grid);

plus_plus_11_assembly = offline_assembly_plus_plus.plus_plus_11;
plus_plus_12_assembly = offline_assembly_plus_plus.plus_plus_12;
plus_plus_21_assembly = offline_assembly_plus_plus.plus_plus_21;
plus_plus_22_assembly = offline_assembly_plus_plus.plus_plus_22;

% for i = 1:1:length(el_subd)
%     para_mapping{i} = poldecomp((para_mapping{i})');
% end

for i = 1:1:length(el_subd)
    para_mapping{i} = eye(2);
end

for i = 1:1:length(el_subd)
    plus_plus_11_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) = ...
        (para_mapping{i}(1,:) * (para_mapping{i}(1,:))') * ...
        plus_plus_11_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) / ...
        det(para_mapping{i})^2;
    plus_plus_12_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) = ...
        (para_mapping{i}(1,:) * (para_mapping{i}(2,:))') * ...
        plus_plus_12_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) / ...
        det(para_mapping{i})^2;
    plus_plus_21_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) = ...
        (para_mapping{i}(2,:) * (para_mapping{i}(1,:))') * ...
        plus_plus_21_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) / ...
        det(para_mapping{i})^2;
    plus_plus_22_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) = ...
        (para_mapping{i}(2,:) * (para_mapping{i}(2,:))') * ...
        plus_plus_22_assembly(gids(el_subd{i},:),gids(el_subd{i},:)) / ...
        det(para_mapping{i})^2;
end

res = plus_plus_11_assembly + plus_plus_12_assembly + ...
    plus_plus_21_assembly + plus_plus_22_assembly;

end