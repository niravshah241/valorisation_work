function [ subd_num ] = get_subdomain_number_extension...
    ( el_num, ref_el_subd )

for i = 1:1:length(ref_el_subd)
    x = find(ref_el_subd{i} == el_num);
    if ~isempty(x)
        break
    end
end

subd_num = i;

end