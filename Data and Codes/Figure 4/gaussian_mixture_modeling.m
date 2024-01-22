clc
clear
%

% cows
load ModelFits/cows_metadata_24_0.975.mat
[nCom nOb] = size(xs_b_train);
options = statset('MaxIter',5000);
for nClst = 1:10
    nClst
    for iter=1:500
        mx = fitgmdist(ZC,nClst,'RegularizationValue',0.000001,'Options',options);
        mdl{nClst,iter} = mx;
        bic(nClst,iter) = mx.BIC;
    end
end
save ModelFits/gaussian_cows_metadata_train_test mdl bic



