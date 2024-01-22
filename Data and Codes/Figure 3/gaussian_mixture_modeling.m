clc
clear
%

% cows
load ModelFits/cows_train_test_16
[nCom nOb] = size(xs_b_train);
thetB = thetB - min(thetB')';
t = ZC*thetB;
[u s v] = svd(t);
ZC = u(:,1:16);thetB = s(1:16,1:16)*v(:,1:16)';
options = statset('MaxIter',5000);

for nClst = 1:10
    nClst
    for iter=1:500
        mx = fitgmdist(ZC,nClst,'RegularizationValue',0.000001,'Options',options);
        mdl{nClst,iter} = mx;
        bic(nClst,iter) = mx.BIC;
    end
end
save ModelFits/gaussian_cows_train_test mdl bic


% % chicken small
% load chicken_small_train_test_10
% [nCom nOb] = size(xs_b_train);
% thetB = thetB - min(thetB')';
% t = ZC*thetB;
% [u s v] = svd(t);
% ZC = u(:,1:10);thetB = s(1:10,1:10)*v(:,1:10)';
% options = statset('MaxIter',5000);
% for nClst = 1:10
%     nClst
%     for iter=1:100
%         mx = fitgmdist(ZC,nClst,'RegularizationValue',0.000001,'Options',options);
%         mdl{nClst,iter} = mx;
%         bic(nClst,iter) = mx.BIC;
%     end
% end
% save gaussian_chicken_small_10 mdl bic


% % chicken
% load chicken_train_test_5
% [nCom nOb] = size(xs_b_train);
% thetB = thetB - min(thetB')';
% t = ZC*thetB;
% [u s v] = svd(t);
% ZC = u(:,1:15);thetB = s(1:15,1:15)*v(:,1:15)';
% options = statset('MaxIter',5000);
% for nClst = 1:10
%     nClst
%     for iter=1:50
%         mx = fitgmdist(ZC,nClst,'RegularizationValue',0.000001,'Options',options);
%         mdl{nClst,iter} = mx;
%         bic(nClst,iter) = mx.BIC;
%     end
% end
% save gaussian_chicken_5 mdl bic


% Human
% load human_train_test_15
% [nCom nOb] = size(xs_b_train);
% thetB = thetB - min(thetB')';
% t = ZC*thetB;
% [u s v] = svd(t);
% ZC = u(:,1:15);thetB = s(1:15,1:15)*v(:,1:15)';
% options = statset('MaxIter',5000);
% 
% for nClst = 1:10
%     nClst
%     for iter=1:50
%         mx = fitgmdist(ZC,nClst,'RegularizationValue',0.000001,'Options',options);
%         mdl{nClst,iter} = mx;
%         bic(nClst,iter) = mx.BIC;
%     end
% end
% save gaussian_human mdl bic
% 
