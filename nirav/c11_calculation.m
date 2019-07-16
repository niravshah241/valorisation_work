c11_ref = c11;
a = 0;
for i = 1:1:length(para_mapping)
    if a < max(max(abs(para_mapping{i})))
        a = max(max(abs(para_mapping{i})));
    end
end
c11 = c11_ref / a;
if (transformed_grid.hmin / grid.hmin) <= a
    disp('correct transformation');
else
    disp('check c11 calculation');
end