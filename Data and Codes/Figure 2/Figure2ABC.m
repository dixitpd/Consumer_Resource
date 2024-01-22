clc
clear
%

subplot(3,3,1)
load ModelFits/cows_157
[nCom nOb] = size(xs_b);
thetB = thetB - min(thetB')';
t = ZC*thetB;
[u s0 v] = svd(t);mx = min(nOb,nCom);
s0 = diag(s0);s0 = s0(1:mx-1);
f = (s0(1:end-1)./s0(2:end));
scatter(1:length(f),f,30,'k','filled')
hold on
plot(1:100,ones(100,1),'k--','linewidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1/2 50])
ylim([5/10 30])

subplot(3,3,4)
hold on
scatter(1:length(s0),s0,30,'k','filled')
t1 = s0(20:80);t2 = log(1:length(t1));
t1 = log(t1);
[a b] = polyfit(t2,t1,1)
xlim([1/2 50])
%ylim([5/10 30])
set(gca,'XScale','log')
set(gca,'YScale','log')

% 

subplot(3,3,2)
clear
load ModelFits/chicken_156
[nCom nOb] = size(xs_b);
thetB = thetB - min(thetB')';
t = ZC*thetB;
[u s0 v] = svd(t);mx = min(nOb,nCom);
s0 = diag(s0);s0 = s0(1:mx-1);
f = (s0(1:end-1)./s0(2:end));
scatter(1:length(f),f,30,'k','filled')
hold on
plot(1:100,ones(100,1),'k--','linewidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1/2 50])
ylim([5/10 30])

subplot(3,3,5)
scatter(1:length(s0),s0,30,'k','filled')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1/2 50])
%ylim([5/10 30])

% 
% 
subplot(3,3,3)
clear
load ModelFits/human_95
[nCom nOb] = size(xs_b);
thetB = thetB - min(thetB')';
t = ZC*thetB;
[u s0 v] = svd(t);mx = min(nOb,nCom);
s0 = diag(s0);s0 = s0(1:mx-1);
f = (s0(1:end-1)./s0(2:end));
scatter(1:length(f),f,30,'k','filled')
hold on
plot(1:100,ones(100,1),'k--','linewidth',1.5)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1/2 50])
ylim([5/10 30])

subplot(3,3,6)
scatter(1:length(s0),s0,30,'k','filled')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1/2 50])
%ylim([5/10 30])