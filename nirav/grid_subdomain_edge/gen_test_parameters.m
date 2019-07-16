function [res] = gen_test_parameters(para_test_1,para_test_2)
para_test_1 = rand(para_test_1(3),1)*(para_test_1(2)-para_test_1(1)) + ...
    para_test_1(1);
para_test_2 = rand(para_test_2(3),1)*(para_test_2(2)-para_test_2(1)) + ...
    para_test_2(1);
[para_test_1,para_test_2] = meshgrid(para_test_1,para_test_2);

para_test_1 = reshape(para_test_1,[],1);
para_test_2 = reshape(para_test_2,[],1);

res = [para_test_1,para_test_2];