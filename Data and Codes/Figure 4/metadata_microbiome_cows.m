clc
clear
%

load ../Data/three_kingdoms_cleaned_up
taxonomy(1,:) = [];
%% Only abundant OTUs in species 1 of cows
sps1 = find(cow_species==1);
ctf_otu      = 0.001;% based on Brian's paper
bacteria_xs  = bacteria_xs(sps1,:);
mn_bact      = mean(bacteria_xs);
goodb        = find(mn_bact > ctf_otu);
restb = 1 - sum(bacteria_xs(:,goodb)')';
xs_b  = [bacteria_xs(:,goodb) restb];
[nSamp nOb] = size(xs_b);
taxonomy_small = taxonomy(goodb,:);
load ../Figure3/ModelFits/test_data_cows.mat

md                       = metadata_cow;
badmets                  = [2 4 5 6 7 8 37];
metnames(badmets)        = [];
md(:,badmets)            = [];
md                       = md(sps1,:);
md_test                  = md(test_data,:);
md_train                 = md;
md_train(test_data,:) = [];

mud      = mean(md_train);
sud      = std(md_train);
md_train = (md_train - mud)./sud;
md_test  = (md_test - mud)./sud;
xs_b_train = xs_b;
xs_b_train(test_data,:) = [];
xs_b_test = xs_b(test_data,:);

nMet = size(md_train,2);
[nSamp nOb] = size(xs_b_train);
% 
% 
nCom = 24;alph   = 0.995;
filen = strcat('ModelFits/cows_metadata_',num2str(nCom),'_',num2str(alph),'.mat')

ZC     = -0.1*rand(nSamp,nCom);
thetB  = 0.1*rand(nCom,nOb);
C      = 0.1*randn(nCom,nMet);

%load(filen)
% hyperparameters
eta = 0.002;
grdnorm = 1;iter = 1;ctf_grad = 0.01;

while grdnorm > ctf_grad
    QB = exp(-[ZC]*thetB);QB = normalize(QB,2,'norm',1);
    deltB = xs_b_train-QB;

    % gradients
    grthetB = [ZC]'*deltB;
    grz     = alph*deltB*thetB' + 2*(1-alph)*(ZC*C - md_train)*C';
    grc     = 2*(1-alph)*ZC'*(ZC*C - md_train);

    % update the variables
    ZC     = ZC - eta*grz;
    thetB  = thetB - eta*grthetB;
    C      = C - eta*grc;

    % errors
    grdnorm = norm(grz)/norm(ZC) + norm(grthetB)/norm(thetB) + norm(grc)/norm(C);

    % Output
    if mod(iter,500) == 0
        save(filen,'taxonomy_small','xs_b_train','ZC','thetB','C','md_test','md_train','QB','xs_b_test','metnames','mud','sud')
        grdnorm
    end
    iter = iter + 1;
end