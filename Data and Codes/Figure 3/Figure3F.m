clc
clear

load ModelFits/gaussian_cows_train_test
load ModelFits/cows_train_test_16
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


load ../Data/three_kingdoms_cleaned_up
sps1         = find(cow_species==1);
ctf_otu      = 0.001;% based on Brian's paper
xs_b_1       = bacteria_xs(sps1,:);
mn_bact      = mean(xs_b_1);
goodb        = find(mn_bact > ctf_otu);
sps2         = find(cow_species==2);
xs_b_2       = bacteria_xs(sps2,goodb);
restb2       = 1 - sum(xs_b_2');
xs_b_2       = [xs_b_2 restb2'];



% nTest = size(xs_b_test,1);
% Ztest = 0.01*randn(nTest,16);
% Z2 = 0.01*randn(size(xs_b_2,1),16);
% iter = 1;etaZ = 0.0001;ctf_grad = 0.01;grdnorm = 1;
% while grdnorm > ctf_grad
%     QB_test = exp(-[Ztest]*thetB);QB_test = normalize(QB_test,2,'norm',1);
%     deltB_test = xs_b_test-QB_test;
%     QB_2 = exp(-[Z2]*thetB);QB_2 = normalize(QB_2,2,'norm',1);
%     deltB_2 = xs_b_2-QB_2;
% 
%     % gradients
%     grz_test     = deltB_test*thetB';
%     grz_2        = deltB_2*thetB';
%     % update the variables
%     Ztest     = Ztest - etaZ*grz_test;
%     Z2        = Z2    - etaZ*grz_2;
%     % errors
%     grdnorm = norm(grz_test)/norm(Ztest)+norm(grz_2)/norm(Z2);
%     % Output
%     if mod(iter,200) == 0
%         grdnorm
%     end
%     iter = iter + 1;
% end
% save('ModelFits/testing_sample_embeddings.mat','Ztest','Z2')
load ModelFits/testing_sample_embeddings.mat

% % 
% % 
% % 
% % %% Embed test data: KL
% subplot(1,2,1)

QB_test = exp(-[Ztest]*thetB);QB_test = normalize(QB_test,2,'norm',1);
QB_2 = exp(-[Z2]*thetB);QB_2 = normalize(QB_2,2,'norm',1);
QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);

shn_test = -nansum( (xs_b_test.*log(QB_test./xs_b_test))' );
shn_2    = -nansum( (xs_b_2.*log(QB_2./xs_b_2))');
shn_trn  = -nansum( (xs_b_train.*log(QB./xs_b_train))' );
p_test   = log10(pdf(gm,Ztest));
p_2      = log10(pdf(gm,Z2));
pTrn     = log10(pdf(gm,ZC));

hold on
scatter(shn_trn,pTrn,25,'k','filled')
scatter(shn_test,p_test,25,'b','filled')
scatter(shn_2,p_2,25,'r','filled')