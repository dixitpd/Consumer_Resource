clc
clear

load ModelFits/gaussian_human_train_test
load ModelFits/human_train_test_16
[nSamp nOb] = size(xs_b_train);
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
% 

%% Species abundance distribution
subplot(1,3,1)
grayColor = [.7 .7 .7];

Zr = random(gm,nSamp);
QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
QB = reshape(QB(:,1:end-1),nSamp*(nOb-1),1);QB(QB < 1e-4) = [];
bx = (-4:0.1:0);
[b a] = hist(log10(QB),bx);b = b/sum(b);b = cumsum(b);b = log10(1-b);
hold on
scatter(a,b,25,'filled','Color',grayColor)
aa = polyfit(a(11:30),b(11:30),1);
x = -4:0.1:0;
y = aa(1)*x+aa(2);
aa(1)
plot(x,y,'Color',grayColor,'linewidth',2.5)

xb = reshape(xs_b_train(:,1:end-1),nSamp*(nOb-1),1);xb(xb<1e-4) = [];
[b a] = hist(log10(xb),bx);b = b/sum(b);b = cumsum(b);b = log10(1-b);
scatter(a,b,25,'k','filled')
aa = polyfit(a(11:30),b(11:30),1);
x = -4:0.1:0;
y = aa(1)*x+aa(2);
plot(x,y,'k--','linewidth',2.5)
ylim([-4 0.5])
aa(1)

% %-1.15: model
% %-1.17: data
% 

%% Shannon
subplot(1,3,2)
clear bx
for iter=1:10
    Zr = random(gm,nSamp);
    QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
    shnQ = -nansum( (QB.*log(QB))' );
    [b a] = ksdensity(shnQ,0:0.025:5);b = b/sum(b);b=b/max(b);
    bx(iter,:) = b;
    ax(iter,:) = a;
end
hold on
mb = mean(bx);sb = std(bx);
ma = mean(ax);
plot(ma,mb+sb,'r')
plot(ma,mb-sb,'r')
inBetween = [mb+sb, fliplr(mb-sb)];
fill([ma fliplr(ma)], inBetween, 'b');
plot(a,mb,'b--')
shnX = -nansum( (xs_b_train.*log(xs_b_train))' );
[b a] = ksdensity(shnX,0:0.025:5);b = b/sum(b);b=b/max(b);
grayColor = [.7 .7 .7];
plot(a,b,'k','linewidth',2.5,'Color',grayColor)
xlim([0 5])
ylim([0 1.1])

% 

%% Bray-Curtis

subplot(1,3,3)
clear a b bx ax
for iter=1:10
    Zr = random(gm,nSamp);
    QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
    fQ = f_braycurtis(QB');fQ(fQ==0) = 1e10;fQ = min(fQ);
    [b a] = ksdensity(fQ,0:0.025:0.5);b = b/sum(b);b=b/max(b);
    bx(iter,:) = b;
    ax(iter,:) = a;
end
hold on
mb = mean(bx);sb = std(bx);
ma = mean(ax);
plot(ma,mb+sb,'b')
plot(ma,mb-sb,'b')
inBetween = [mb+sb, fliplr(mb-sb)];
fill([ma fliplr(ma)], inBetween, 'b');

fX = f_braycurtis(xs_b_train');fX(fX==0) = 1e10;fX = min(fX);
[b a] = ksdensity(fX,0:0.025:0.5);b = b/sum(b);b=b/max(b);
plot(a,b,'k','linewidth',2.5,'Color',grayColor)
ylim([0 1.1])
% 
clear ax bx
for iter=1:10
    Zr = random(gm,nSamp);
    QB = exp(-Zr*thetB);QB = normalize(QB,2,'norm',1);
    fQ = f_braycurtis(QB');fQ(fQ==0) = [];
    [b a] = ksdensity(fQ,0:0.025:2);b = b/sum(b);b=b/max(b);
    bx(iter,:) = b;
    ax(iter,:) = a;
end
mb = mean(bx);sb = std(bx);ma = mean(ax);
ma = 0.5+ma;
plot(ma,mb+sb,'b')
plot(ma,mb-sb,'b')
inBetween = [mb+sb, fliplr(mb-sb)];
fill([ma fliplr(ma)], inBetween, 'b');
fX = f_braycurtis(xs_b_train');fX(fX==0) = [];
[b a] = ksdensity(fX,0:0.025:2);b = b/sum(b);b=b/max(b);
a = a + 0.5;
plot(a,b,'k','linewidth',2.5,'Color',grayColor)
xlim([0 2])


