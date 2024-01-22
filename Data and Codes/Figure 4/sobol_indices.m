clc
clear


load ModelFits/gaussian_cows_metadata_train_test
load ModelFits/cows_metadata_24_0.975
[nSamp nOb] = size(xs_b_train);
nMet = size(md_train,2);

mbic = mean(bic');
ii = find(mbic==min((mbic)));
jj = find(min(bic(ii,:))==bic(ii,:));
mdlx = mdl{ii,jj};
gm = gmdistribution(mdlx.mu,mdlx.Sigma)


for metid = 1:nMet
    metid

    clear mQ pH xa QM Zr
    mi = -3;
    mx = 3;
    nC = 40;
    dx = (mx-mi)/nC;


    Zr = random(gm,100000);
    QM = Zr*C;strx = QM(:,metid);
    clear mQ xa mb
    for xx = 1:25
        clear bc
        for i=1:nC
            st = mi + dx*(i-1);
            en = st + dx;
            xa(i) = 0.5*(st+en);
            ii = find((strx > st).*(strx < en));
            if length(ii) > 100
                ii = randsample(ii,100);
            end
            if length(ii) > 10
                QB = exp(-Zr(ii,:)*thetB);QB = normalize(QB,2,'norm',1);
                g = f_braycurtis(QB');g(g==0) = [];
                bc(i) = mean(g);
            else
                bc(i) = nan;
            end
        end
        mb(xx,:) = bc;
    end
    mbc(metid,:) = mean(mb);
    sbc(metid,:) = std(mb);

end
save('ModelFits/sobol_behavior_cows','mbc','sbc','xa')

