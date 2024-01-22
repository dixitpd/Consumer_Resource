clc
clear

load ../Figure4/ModelFits/gaussian_cows_metadata_train_test
load ../Figure4/ModelFits/cows_metadata_24_0.975.mat
[nSamp nOb] = size(xs_b_train);
nMet = size(md_train,2);

mbic = mean(bic');
ii = find(mbic==min((mbic)));
jj = find(min(bic(ii,:))==bic(ii,:));
mdlx = mdl{ii,jj};
gm = gmdistribution(mdlx.mu,mdlx.Sigma);
load ../Figure4/ModelFits/sobol_behavior_cows
for xx=1:50
    Zr = random(gm,nSamp);
    QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
    g = f_braycurtis(QB');x = g;x(x==0) = 1e10;x = min(x);
    g(g==0) = [];
    pop_var(xx) = mean(g);
    pop_min(xx) = mean(x);
end
mpop_var = mean(pop_var);
spop_var = std(pop_var);
mpop_min = mean(pop_min);
spop_min = std(pop_min);

nC = length(xa);

[a b] = sort(nanmean( ((mbc/mpop_var - 1))' ) );bx = b;