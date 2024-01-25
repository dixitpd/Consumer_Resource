clc
clear

load ../Figure4/ModelFits/gaussian_cows_metadata_train_test
load ../Figure4/ModelFits/cows_metadata_24_0.975
[nSamp nOb] = size(xs_b_train);
nMet = size(md_train,2);
nK = size(ZC,2);

mbic = mean(bic');
ii = find(mbic==min((mbic)));
jj = find(min(bic(ii,:))==bic(ii,:));
mdlx = mdl{ii,jj};
gm = gmdistribution(mdlx.mu,mdlx.Sigma)

id1 = 6;id2 = 23; 
ct  = 1;nIter = 250000;
gamm = 0;
Z0 = [];
for iz = 1:100
    acp = 0;
    ZO  = random(gm,1);ZN = ZO;eO = 1e100;
    for iter=1:nIter
        ZN = ZO + 0.05*randn(1,nK);
        M  = ZN*C;m1 = M(id1);m2 = M(id2);
        eN = -log(pdf(gm,ZN)) + gamm*(m1-3)*(m1-3) + gamm*(m2-3)*(m2-3);
        if exp(-(eN-eO)) > rand
            eO = eN;
            ZO = ZN;
            acp = acp + 1;
        end

        if iter > 0.5*nIter && mod(iter,100) == 0
            M = ZO*C;
            mx(ct,1) = M(id1);mx(ct,2) = M(id2);
            Z0(ct,:) = ZO;
            ct = ct + 1;
        end
    end
    [iz acp/nIter]
end
QB0 = exp(-Z0*thetB);QB0 = normalize(QB0,2,'norm',1);

scatter(md_train(:,id1),md_train(:,id2),50,'k','filled')
hold on
scatter(mx(:,1),mx(:,2),50,'k')

ct  = 1;
gamm = 2;
Zx = [];
for iz = 1:100
    acp = 0;
    ZO  = random(gm,1);ZN = ZO;eO = 1e100;
    for iter=1:nIter
        ZN = ZO + 0.05*randn(1,nK);
        M  = ZN*C;m1 = M(id1);m2 = M(id2);
        eN = -log(pdf(gm,ZN)) + gamm*(m1-3)*(m1-3) + gamm*(m2-3)*(m2-3);
        if exp(-(eN-eO)) > rand
            eO = eN;
            ZO = ZN;
            acp = acp + 1;
        end

        if iter > 0.5*nIter && mod(iter,100) == 0
            M = ZO*C;
            gx(ct,1) = M(id1);gx(ct,2) = M(id2);
            Zx(ct,:) = ZO;
            ct = ct + 1;
        end
    end
    [iz acp/nIter]
end
scatter(gx(:,1),gx(:,2),50,'b','filled')
QB = exp(-Zx*thetB);QB = normalize(QB,2,'norm',1);

% [65   114   103   148   108]
% [65   114   103    31    69]