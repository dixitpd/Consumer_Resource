clc
clear
% 

load ModelFits/gaussian_cows_metadata_train_test
load ModelFits/cows_metadata_24_0.975.mat
[nSamp nOb] = size(xs_b_train);
nMet = size(md_train,2);

mbic = mean(bic');
ii = find(mbic==min((mbic)));
jj = find(min(bic(ii,:))==bic(ii,:));
mdlx = mdl{ii,jj};
gm = gmdistribution(mdlx.mu,mdlx.Sigma);
load ModelFits/sobol_behavior_cows
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
b1 = [b(1:5)]
list_cows = b;
save('ModelFits/list_cows','list_cows')
% 
% 
for metid = 1:3
    subplot(1,3,metid)
    hold on
    ylim([0 1.25])
    plot(xa,ones(nC,1),'k--')
    errorbar(xa*sud(b(metid)) + mud(b(metid)),mbc(b(metid),:)/mpop_var,sbc(b(metid),:)/mpop_var,'b','linewidth',1.5)
    fontsize(15,"points")
end

figure

%
load ModelFits/cows_fits
% ranksum(bc_good,bc_all)
% ranksum(bc_good,ave_bc_data)
% ranksum(ave_bc_data,bc_all)

subplot(1,4,1)
hold on
ylim([0 1])
boxplot(bc_good)
scatter(ones(size(bc_good)).*(1+(rand(size(bc_good))-0.5)/10),bc_good,5,'k','filled')
xlim([0.8 1.2])

subplot(1,4,2)
hold on
ylim([0 1])
boxplot(bc_bot)
scatter(ones(size(bc_bot)).*(1+(rand(size(bc_bot))-0.5)/10),bc_bot,5,'k','filled')
xlim([0.8 1.2])

subplot(1,4,3)
hold on
ylim([0 1])
boxplot(bc_all)
scatter(ones(size(bc_all)).*(1+(rand(size(bc_all))-0.5)/10),bc_all,5,'k','filled')
xlim([0.8 1.2])

subplot(1,4,4)
hold on
ylim([0 1])
boxplot(ave_bc_data)
ave_bc_data = randsample(ave_bc_data,100);
scatter(ones(size(ave_bc_data)).*(1+(rand(size(ave_bc_data))-0.5)/10),ave_bc_data,5,'k','filled')
xlim([0.8 1.2])



