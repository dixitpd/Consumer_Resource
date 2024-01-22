clc
clear

load ../Figure2/ModelFits/cows_16
[nSamp nOb] = size(xs_b);
thetB = thetB - min(thetB')';
t = ZC*thetB;
[u s v] = svd(t);
ZC = u(:,1:16);thetB = s(1:16,1:16)*v(:,1:16)';


load ../Data/three_kingdoms_cleaned_up
taxonomy(1,:) = [];
sps1 = find(cow_species==1);
md                       = metadata_cow;
badmets                  = [2 4 5 6 7 8 37];
metnames(badmets)        = [];
md(:,badmets)            = [];
md                       = md(sps1,:);

for iter=1:size(md,2)
    [a b] = lasso(ZC,md(:,iter),'lambda',0);
    intr  = b.Intercept;
    [crx(iter) prx(iter)] = corr(ZC*a + intr,md(:,iter));
end
