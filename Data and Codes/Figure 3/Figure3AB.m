clc
clear

load ModelFits/gaussian_cows
load ../Figure2/ModelFits/cows_16
[nSamp nOb] = size(xs_b);
thetB = thetB - min(thetB')';
t = ZC*thetB;
[u s v] = svd(t);
ZC = u(:,1:16);thetB = s(1:16,1:16)*v(:,1:16)';

mbic = mean(bic');
ii = find(mbic == min(mbic));
jj = find(bic(ii,:)==min(bic(ii,:)));
mdlx = mdl{ii,jj};
gm = gmdistribution(mdlx.mu,mdlx.Sigma)
% 

%% Means and variances

subplot(1,3,1)
Zr = random(gm,nSamp);
qb = exp(-Zr*thetB);qb = normalize(qb,2,'norm',1);
for iter=1:size(qb,1)
    QB(iter,:) = mnrnd(5000,qb(iter,:));
end
QB = normalize(QB,2,'norm',1);
mQ = mean(QB(:,1:end-1));sQ = std(QB(:,1:end-1));
mX = mean(xs_b(:,1:end-1));sX = std(xs_b(:,1:end-1));
hold on
scatter(log10(mX),log10(mQ),50,'k','filled')
scatter(log10(sX),log10(sQ),50,'ko','linewidth',1.5)
plot([-3.5 -0.5 ],[-3.5 -0.5],'r--','linewidth',1.5)
xlim([[-3.5 -0.5 ]])
ylim([-3.5 -0.5 ])

% 

%% Covariance

subplot(1,3,2)
Zr = random(gm,nSamp);
QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
hold on
plot([1e-10 1],[1e-10 1],'r--')
cQ = cov(QB(:,1:end-1));cX = cov(xs_b(:,1:end-1));cQ = cQ - diag(diag(cQ));cX = cX - diag(diag(cX));
cQ(cQ==0) = [];cX(cX==0) = [];
scatter(log10(abs(cX)),log10(abs(cQ)),30,'filled','k')
plot([-10 -2],[-10 -2],'r--','linewidth',1.5)
ylim([-10 -2])
xlim([-10 -2])
% 
% 

subplot(1,3,3)
x = [1 2 3 4];
my = [0.9987 0.9349 0.4484 0.5691]
sy = [5.7285e-04 0.0122 0.1217 0.1813]
hold on
bar(x,my,'b')
errorbar(x,my,sy,'ko','linewidth',1.5)

% for iter=1:10
%     Zr = random(gm,nSamp);
%     QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
%     mQ = mean(QB(:,1:end-1));mX = mean(xs_b(:,1:end-1));
%     cc(iter) = corr((mQ)',(mX)');
% end

% for iter=1:10
%     Zr = random(gm,nSamp);
%     QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
%     cQ = cov(QB(:,1:end-1));cX = cov(xs_b(:,1:end-1));cQ = cQ - diag(diag(cQ));cX = cX - diag(diag(cX));
%     cQ(cQ==0) = [];cX(cX==0) = [];
%     cc(iter) = corr((cX)',(cQ)');
% end

% 
% %% 3-body correlations
% for iter=1:10
%     iter
%     Zr = random(gm,nSamp);
%     QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
% 
%     for ix=1:100000
%         ijk = randsample(nOb-1,3);
%         i = ijk(1);j = ijk(2);k = ijk(3);
%         c3Q(ix) = mean((QB(:,i)-mQ(i)).*(QB(:,j)-mQ(j)).*(QB(:,k)-mQ(k)));
%         c3X(ix) = mean((xs_b(:,i)-mX(i)).*(xs_b(:,j)-mX(j)).*(xs_b(:,k)-mX(k)));
%     end
%     cc(iter) = corr(c3Q',c3X');
% 
% end
% % 
%% 4-body correlations

% for iter=1:10
%     iter
%     Zr = random(gm,nSamp);
%     QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
% 
%     for ix=1:50000
%         ijkl = randsample(nOb-1,4);
%         i = ijkl(1);j = ijkl(2);k = ijkl(3);l = ijkl(4);
%         c4Q(ix) = mean((QB(:,i)-mQ(i)).*(QB(:,j)-mQ(j)).*(QB(:,k)-mQ(k)).*(QB(:,l)-mQ(l)));
%         c4X(ix) = mean((xs_b(:,i)-mX(i)).*(xs_b(:,j)-mX(j)).*(xs_b(:,k)-mX(k)).*(xs_b(:,l)-mX(l)));
%     end
%     cc(iter) = corr(c4Q',c4X');
% 
% end
% 

