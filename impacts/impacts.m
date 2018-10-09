% Load Earth Impact Database

eid = importdataset('EID.csv','\t');
eid.Age = (eid.Age_Max+eid.Age_Min)/2;
eid.Age_Sigma = (eid.Age_Max-eid.Age_Min)/4;


%% Load crustal area distribution
load GSCmap.mat % km^2 bedrock area per Myr

figure; plot(GSCmap.Age,GSCmap.Area,'LineWidth',2)
set(gca,'xdir','reverse')
xlabel('Age (Ma)');
set(gca,'yscale','log')
ylabel('Exposed Continental Surface Area (km^2 / yr)')
ylim([5*10^3,10^7])
formatfigure
warning('off', 'MATLAB:print:FigureTooLargeForPage')
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig, 'ContinentalSurfaceArea.pdf','pdf')

% Plot area distribution
figure; plot(GSCmap.Age,cumsum(GSCmap.Area,1,'reverse')./sum(GSCmap.Area),'LineWidth',1)
xlabel('Age (Ma)'); ylabel('Cumulative continental area (fraction of total)')
set(gca,'Xdir','Reverse')
formatfigure;
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig,'CumulativeAreaFraction.pdf','pdf')


%% Plot terrestrial craters vs exposed surface area on the same plot
binedges = 0.01:100:4000.01;
bincenters = center(binedges);
c = lines(7);

figure;
yyaxis('left')

% t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>2;
% [Ne,~] = histcounts(eid.Age(t),binedges);
% bar(bincenters,Ne,1,'FaceColor',c(5,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>10;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne,1,'FaceColor',c(6,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>20;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne,1,'FaceColor',c(1,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>100;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne,1,'FaceColor',c(4,:))

ylabel('Impact Craters / 10^8 yr')
ylim([0,25.5]);
set(gca,'ytick',0:5:25)

hold on;
plot([717.4,717.4],get(gca,'ylim'),'Color',c(1,:),'LineWidth',1) %Sturtian

yyaxis('right')
plot(GSCmap.Age,GSCmap.Area,'LineWidth',2)
set(gca,'xdir','reverse')
xlabel('Age (Ma)');
ylabel('Exposed Area (km^2 / Myr)')
ylim([0,5*10^5]);
set(gca,'ytick',(0:5)*1E5)

xlim([0 2500])
formatfigure
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig, 'TerrestrialCratersvsSurfaceArea.pdf','pdf')

%% Plot terrestrial craters per exposed surface area

binedges = 0.01:100:4000.01;
bincenters = center(binedges);

CrustArea = binsum(GSCmap.Age,GSCmap.Area,binedges);
CrustArea = CrustArea/10^7; % Normalize per 10^7 km^2
figure;
c = lines(7);

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>2;
[Ne,~] = histcounts(eid.Age(t),binedges);
figure; bar(bincenters,Ne./CrustArea,1,'FaceColor',c(5,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>10;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne./CrustArea,1,'FaceColor',c(6,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>20;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne./CrustArea,1,'FaceColor',c(1,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>100;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne./CrustArea,1,'FaceColor',c(4,:))

xlabel('Age (Ma)')
ylabel('Craters per 100 Myr per 10^7km^2')
set(gca,'XDir','reverse')
% set(gca,'ytick',0:20:100)
formatfigure;

hold on;
plot([717.4,717.4],get(gca,'ylim'),'Color',c(1,:),'LineWidth',1) %Sturtian


legend({'>2 km','>10 km','>20 km','>100 km'},'Location','Northwest')
xlim([0 1600])
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig, 'TerrestrialCratersPerExposedArea.pdf','pdf')


%% Plot terrestrial craters per cumulative surface area

CrustAreaCumulative = cumsum(CrustArea,2,'reverse');
CrustAreaCumulative = CrustAreaCumulative/10^7; % Normalize per 10^7 km^2

figure;
c = lines(7);

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>2;
[Ne,~] = histcounts(eid.Age(t),binedges);
figure; bar(bincenters,Ne./CrustAreaCumulative,1,'FaceColor',c(5,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>10;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne./CrustAreaCumulative,1,'FaceColor',c(6,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>20;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne./CrustAreaCumulative,1,'FaceColor',c(1,:))

t = eid.Age>0 & eid.Age_Sigma<75 & eid.Diameter>100;
[Ne,~] = histcounts(eid.Age(t),binedges);
hold on; bar(bincenters,Ne./CrustAreaCumulative,1,'FaceColor',c(4,:))

xlabel('Age (Ma)')
ylabel('Craters / 100 Myr / 10^7km^2 cumulative area')
set(gca,'XDir','reverse')
% set(gca,'ytick',0:20:100)
formatfigure;

hold on;
plot([717.4,717.4],get(gca,'ylim'),'Color',c(1,:),'LineWidth',1) %Sturtian


legend({'>2 km','>10 km','>20 km','>100 km'},'Location','Northwest')
xlim([0 1600])
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig, 'TerrestrialCratersPerCumulativeArea.pdf','pdf')


