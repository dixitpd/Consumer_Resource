clc
clear
%
load ../Figure4/ModelFits/gaussian_cows_metadata_train_test
load ../Figure4/ModelFits/cows_metadata_24_0.975
[nSamp nOb] = size(xs_b_train);
nMet = size(md_train,2);
nK = size(ZC,2);

subplot(2,2,1)
load ModelFits/propionate
Z0 = Z0(1:250:end,:);Zx = Zx(1:250:end,:);
M0 = Z0*C;M0 = M0.*sud + mud;
Mx = Zx*C;Mx = Mx.*sud + mud;
Md = md_train.*sud + mud;
id1 = 26; % Starch
id2 = 23; % pH

gridx1 = 150:10:400;l1 = length(gridx1);
gridx2 = 5.5:0.1:8;l2 = length(gridx2);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];

x = [M0(:,id1) M0(:,id2)];
[b a] = ksdensity(x,xi);
xa = reshape(a(:,1),l2,l1);
ya = reshape(a(:,2),l2,l1);
b = reshape(b,l2,l1);b(b<1e-3) = nan;
[m n] = contour(xa,ya,b,'k--');
n.LineWidth = 1;

hold on

x = [Md(:,id1) Md(:,id2)];
[b a] = ksdensity(x,xi);
xa = reshape(a(:,1),l2,l1);
ya = reshape(a(:,2),l2,l1);
b = reshape(b,l2,l1);b(b<1e-3) = nan;
[m n] = contour(xa,ya,b,'k');
n.LineWidth = 1;

x = [Mx(:,id1) Mx(:,id2)];
[b a] = ksdensity(x,xi);
xa = reshape(a(:,1),l2,l1);
ya = reshape(a(:,2),l2,l1);
b = reshape(b,l2,l1);b(b<1e-3) = nan;
[m n] = contour(xa,ya,b,'r');
n.LineWidth = 1;
ylim([5 8])
% 
subplot(2,2,2)
hold on
load ModelFits/propionate
M0 = Z0*C;M0 = M0.*sud + mud;
Mx = Zx*C;Mx = Mx.*sud + mud;
Md = md_train.*sud + mud;
m0 = mean(M0(:,25:30));
s0 = std(M0(:,25:30));
mx = mean(Mx(:,25:30));
sx = mean(Mx(:,25:30));
bar(3*(1:6)+1.1,m0,0.25,'k')
bar(3*(1:6)+1.9,mx,0.25,'r')
hold on


% subplot(2,3,3)
% hold on
% load ModelFits/propionate
% QB = exp(-Zx*thetB);QB = normalize(QB,2,'norm',1);
% QB0 = exp(-Z0*thetB);QB0 = normalize(QB0,2,'norm',1);
% m0 = mean(QB0(:,1:end-1));mx = mean(QB(:,1:end-1));
% scatter(m0,mx,20,'k','filled')
% xlim([1e-4 1])
% ylim([1e-6 1])
% plot([1e-4 1],[1e-4 1],'r--','LineWidth',2)
% plot([1e-4 1],[2e-4 2],'k--')
% plot([2e-4 2],[1e-4 1],'k--')
% plot([1e-4 1],[1e-2 1e-2],'k--')
% plot([1e-2 1e-2],[1e-4 1],'k--')
% set(gca,'xScale','log')
% set(gca,'yScale','log')




subplot(2,2,3)
load ModelFits/sara
Z0 = Z0(1:250:end,:);Zx = Zx(1:250:end,:);
M0 = Z0*C;
Mx = Zx*C;
Md = md_train;
id1 = 6; % Starch
id2 = 23; % pH

gridx1 = 1:0.1:10;l1 = length(gridx1);
gridx2 = 5.5:0.1:8;l2 = length(gridx2);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];

x = [M0(:,id1)*sud(id1) + mud(id1) M0(:,id2)*sud(id2) + mud(id2)];
[b a] = ksdensity(x,xi);
xa = reshape(a(:,1),l2,l1);
ya = reshape(a(:,2),l2,l1);
b = reshape(b,l2,l1);b(b<1e-3) = nan;
[m n] = contour(xa,ya,b,'k--');
n.LineWidth = 1;

hold on

x = [Md(:,id1)*sud(id1) + mud(id1) Md(:,id2)*sud(id2) + mud(id2)];
[b a] = ksdensity(x,xi);
xa = reshape(a(:,1),l2,l1);
ya = reshape(a(:,2),l2,l1);
b = reshape(b,l2,l1);b(b<1e-3) = nan;
[m n] = contour(xa,ya,b,'k');
n.LineWidth = 1;

x = [Mx(:,id1)*sud(id1) + mud(id1) Mx(:,id2)*sud(id2) + mud(id2)];
[b a] = ksdensity(x,xi);
xa = reshape(a(:,1),l2,l1);
ya = reshape(a(:,2),l2,l1);
b = reshape(b,l2,l1);b(b<1e-3) = nan;
[m n] = contour(xa,ya,b,'r');
n.LineWidth = 1;
xlim([0 10])
ylim([5 8])


subplot(2,2,4)
hold on
load ModelFits/sara
QB = exp(-Zx*thetB);QB = normalize(QB,2,'norm',1);
QB0 = exp(-Z0*thetB);QB0 = normalize(QB0,2,'norm',1);
m0 = mean(QB0(:,1:end-1));mx = mean(QB(:,1:end-1));
scatter(m0,mx,20,'k','filled')
xlim([1e-4 1])
ylim([1e-4 1])
plot([1e-4 1],[1e-4 1],'r--','LineWidth',2)
plot([1e-4 1],[2e-4 2],'k--')
plot([2e-4 2],[1e-4 1],'k--')
plot([1e-4 1],[1e-2 1e-2],'k--')
plot([1e-2 1e-2],[1e-4 1],'k--')
set(gca,'xScale','log')
set(gca,'yScale','log')

