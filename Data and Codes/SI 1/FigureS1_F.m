clc
clear

load ModelFits/gaussian_chicken_train_test
load ModelFits/chicken_train_test_16
[nSamp nOb] = size(xs_b_train);
thetB = thetB - min(thetB')';
t = ZC*thetB;
[u s v] = svd(t);
ZC = u(:,1:16);thetB = s(1:16,1:16)*v(:,1:16)';
% 

mbic = mean(bic');
ii = find(mbic == min(mbic));
jj = find(bic(ii,:)==min(bic(ii,:)));
mdlx = mdl{ii,jj};
gm = gmdistribution(mdlx.mu,mdlx.Sigma)



% nTest = size(xs_b_test,1);
% Ztest = 0.01*randn(nTest,16);
% iter = 1;etaZ = 0.0001;ctf_grad = 0.01;grdnorm = 1;
% while grdnorm > ctf_grad
%     QB_test = exp(-[Ztest]*thetB);QB_test = normalize(QB_test,2,'norm',1);
%     deltB_test = xs_b_test-QB_test;
% 
%     % gradients
%     grz_test     = deltB_test*thetB';
%     % update the variables
%     Ztest     = Ztest - etaZ*grz_test;
%     % errors
%     grdnorm = norm(grz_test)/norm(Ztest);
%     % Output
%     if mod(iter,200) == 0
%         grdnorm
%     end
%     iter = iter + 1;
% end
% save('ModelFits/testing_sample_embeddings.mat','Ztest')
load ModelFits/testing_sample_embeddings.mat

% % 
% % 
% % 
% % %% Embed test data: KL
% subplot(1,2,1)

QB_test = exp(-[Ztest]*thetB);QB_test = normalize(QB_test,2,'norm',1);
QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);

shn_test = -nansum( (xs_b_test.*log(QB_test./xs_b_test))' );
shn_trn  = -nansum( (xs_b_train.*log(QB./xs_b_train))' );
p_test   = log10(pdf(gm,Ztest));
pTrn     = log10(pdf(gm,ZC));

hold on
scatter(shn_trn,pTrn,25,'k','filled')
scatter(shn_test,p_test,25,'b','filled')
