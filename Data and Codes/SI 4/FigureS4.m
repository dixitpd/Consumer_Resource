clc
clear
% 

load ../Figure4/ModelFits/gaussian_cows_metadata_train_test
load ../Figure4/ModelFits/cows_metadata_24_0.975.mat
[nSamp nOb] = size(xs_b_train);
nMet = size(md_train,2);

mbic = mean(bic');
ii = find(mbic==min((mbic)));
jj = find(min(bic(ii,:))==bic(ii,:));
mdlx = mdl{ii,jj};
gm = gmdistribution(mdlx.mu,mdlx.Sigma);

Zr = random(gm,nSamp);
QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
M  = Zr*C;

t1 = corr(QB,M);t1 = reshape(t1,nMet*nOb,1);
t2 = corr(xs_b_train,md_train);t2 = reshape(t2,nMet*nOb,1);
t3 = corr(xs_b_test,md_test);t3 = reshape(t3,nMet*nOb,1);

hold on
scatter(t2,t3,15,'k','filled')
scatter(t2,t1,15,'r','filled')
plot([-1 1],[-1 1],'k--')

