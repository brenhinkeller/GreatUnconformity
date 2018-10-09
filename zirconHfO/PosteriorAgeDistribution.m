load zircon
load mczircon

%%

[a,b]=hist(zircon.Age(zircon.Age>0&zircon.Age<4500),100);
figure; bar(b,a,1)
hold on; plot(linspace(0,4500,100),ones(1,100).*max(a),'k','LineWidth',1.5)
set(gca,'YtickLabel','')
ylim([0 1.1*max(a)])
xlabel('Age (Ma)')
ylabel('Number of samples')
title('Prior age distribution')
formatfigure
saveas(gcf,'Prior age distribution','pdf')

[a,b]=hist(mczircon.Age(mczircon.Age>0&mczircon.Age<4500),100);
figure; bar(b,a,1)
hold on; plot(linspace(0,4500,100),ones(1,100).*max(a),'k','LineWidth',1.5)
set(gca,'YtickLabel','')
ylim([0 1.1*max(a)])
xlabel('Age (Ma)')
ylabel('Number of samples')
title('Posterior age distribution')
formatfigure
saveas(gcf,'Posterior age distribution','pdf')
