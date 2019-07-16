function [ subd_num ] = get_subdomain_number( el_num, ref_el_subd )
%GET_SUBDOMAIN_NUMBER Summary of this function goes here
%   Detailed explanation goes here

for i = 1:1:length(ref_el_subd)
    x = find(ref_el_subd{i} == el_num);
    if ~isempty(x)
        break
    end
end

subd_num = i;

end