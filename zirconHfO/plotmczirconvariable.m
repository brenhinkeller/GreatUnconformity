%% Single element, just zircon

if ~exist('mczircon','var')
    load mczircon
end
if ~exist('zircon','var')
    load zircon
end

Elem='eHf_initial';
agemin=0;
agemax=4350;
nbins = 145;
binoverlap = 3;
% nbins = 124;
% binoverlap = 0;

test=~isnan(mczircon.(Elem));
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
figure; errorbar(c,m,2*e,'.')
xlabel('Age (Ma)'); ylabel(Elem); xlim([agemin, agemax])
set(gca,'xdir','reverse');
xlim([0 4500]);
ylim([-10 10]);
set(gca,'ytick',[-10 -5 0 5 10])
formatfigure

agebincenters=c;
ehf=m;

%%
Elem='Lu176_Hf177';
agemin=0;
agemax=4200;
nbins = 140;
test=~isnan(mczircon.(Elem));
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,3);
figure; errorbar(c,m,2*e,'.r')
xlabel('Age (Ma)'); ylabel(Elem); xlim([agemin, agemax])
set(gca,'xdir','reverse');
xlim([0 4500]);
formatfigure


%%
Elem='d18O';
agemin=0;
agemax=4360;
nbins = 109;
test=~isnan(mczircon.(Elem));
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,3);
figure; errorbar(c,m,2*e,'.r')
xlabel('Age (Ma)'); ylabel(Elem); xlim([agemin, agemax])
set(gca,'xdir','reverse');
xlim([0 4500]);
ylim([4 10])
formatfigure

%% Cumulative
Elem='d18O';
agemin=0;
agemax=4360;
nbins = 109;
test=~isnan(mczircon.(Elem));
[c,m,e]=bincumulative(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins);
figure; errorbar(c,m,2*e,'.r')
xlabel('Age (Ma)'); ylabel(['Cumulative igneous ' Elem]); xlim([agemin, agemax])
set(gca,'xdir','reverse');
xlim([0 4500]);
ylim([5 7])
formatfigure

% Prediction for hydrosphere

mHydrosphere = 1386E6*(10^3)^3*1000;
mCrust = 3.77E22; % From CRUST 2.0 (including sediments)
mIgneousCrust = 3.77E22*(1-0.016-0.038); % From CRUST 2.0

d18OMantle = 5.18; % +/- 0.28, 2-sigma (Mattey et al, 1994)
d18OMantleSigma = 0; %0.14;

d18OComplementary = d18OMantle + (d18OMantle - m)*mIgneousCrust/mHydrosphere*0.45./(16/18);

figure; errorbar(c,d18OComplementary,2*sqrt(e.^2+d18OMantleSigma.^2)*mIgneousCrust/mHydrosphere*0.45./(16/18),'.r')
set(gca,'xdir','reverse');
xlabel('Age (Ma)'); ylabel('Complementary hydrosphere d18O');
xlim([0 4500]);
formatfigure



%% Resample and plot age spectrum (no age weighting)

err=zircon.err2srel.Age;
err(isnan(err))=nanmean(err); % Assign average value to missing uncerts
test=err./zircon.Age < 50;
err(test) = 50./zircon.Age(test); % Set absolute minimum uncert
ages=bsresample(zircon.Age,err,10^7);

figure; hist(ages,1000)


%% Crust subduction signature

agemin=0;
agemax=4350;
nbins = 145;

Elem='eHf_initial';
test=~isnan(mczircon.(Elem));
[c,m1,e1]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,3);

Elem='d18O';
test=~isnan(mczircon.(Elem));
[c,m2,e2]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,3);

m1 = detrend(m1);
m2 = detrend(m2);

e1 = e1./nanstd(m1);
m1 = (-m1+nanmean(m1))./nanstd(m1);
e2 = e2./nanstd(m2);
m2 = (m2-nanmean(m2))./nanstd(m2);

figure; errorbar(c,m1,2*e1,'.'); xlim([0 4500]); xlabel('Age (Ma)'); ylabel('eHf signal, standardized'); set(gca,'xdir','reverse'); formatfigure;
figure; errorbar(c,m2,2*e2,'.'); xlim([0 4500]); xlabel('Age (Ma)'); ylabel('d18o signal, standardized'); set(gca,'xdir','reverse'); formatfigure;


m = nanmean([m1;m2],1);
e = nanmean([e1;e2],1);

e = e./nanstd(m);
m = (m-min(m))./nanstd(m);
figure; errorbar(c,m,2*e,'.r')

xlabel('Age (Ma)'); ylabel('Sediment subduction signature'); xlim([agemin agemax])
set(gca,'xdir','reverse')
formatfigure;


%%
figure; plot(c(4:end-3),m(1:end-6)-m(7:end))
xlabel('Age (Ma)'); ylabel('Sediment subduction differential'); xlim([agemin agemax])
set(gca,'xdir','reverse')
formatfigure;

%% Crust subduction covariance

agemin=0;
agemax=4350;
% nbins = 145;
% binoverlap = 3;
nbins = 290;
binoverlap = 6;

% Calculate zircon eHf timeseries
Elem='eHf_initial';
test=~isnan(mczircon.(Elem));
[c,ehf,ehf_sigma]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);

% Calculate ziron eHf timeseries
Elem='d18O';
test=~isnan(mczircon.(Elem));
[c,d18o,d18o_sigma]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);

% De-trend and unit variance
d18o = detrend(d18o);
d18o = d18o/std(d18o);
ehf = detrend(ehf);
ehf = -ehf/std(ehf); % Negative since lower eHf means more crust remelting
% figure; plot(c,ehf)
% figure; plot(c,d18o)

% Average signal
sigaverage = nanmean([ehf; d18o],1);
sigaverage = sigaverage/std(sigaverage);
figure; plot(c,sigaverage)
xlabel('Age (Ma)'); ylabel('Average Signal')
set(gca,'xdir','reverse')
formatfigure

% Find times of high covariance between eHf and d18o
windowwidth = 8;
cs = c((1+windowwidth/2):(end-windowwidth/2)); % Bin centers for moving window
covarianceHfO = NaN(1,nbins-windowwidth); % covariance
sehf = NaN(1,nbins-windowwidth);
sd18o = NaN(1,nbins-windowwidth);
saverage = NaN(1,nbins-windowwidth);
for i=1:(nbins-windowwidth)
    % Calculate and store the covariance
    temp = cov(ehf(i:i+windowwidth),d18o(i:i+windowwidth));
    covarianceHfO(i) = temp(1,2);
    
    % Fit a line through the window to determine slope and intercept
    pehf = polyfit(c(i:i+windowwidth),ehf(i:i+windowwidth),1);
    pd18o = polyfit(c(i:i+windowwidth),d18o(i:i+windowwidth),1);
    paverage = polyfit(c(i:i+windowwidth),sigaverage(i:i+windowwidth),1);
    
    % Extract and store the slope (minus sign is because time flows towards zero Ma)
    sehf(i) = -pehf(1);
    sd18o(i) = -pd18o(1);
    saverage(i) = -paverage(1);
end

% figure; plot(cs, sehf); xlabel('Age (Ma)'); ylabel('eHf signal, standardized'); set(gca,'ydir','reverse','xdir','reverse'); formatfigure;
% figure; plot(cs, sd18o); xlabel('Age (Ma)'); ylabel('d18o signal, standardized'); set(gca,'xdir','reverse'); formatfigure;

% Average slope
aves = nanmean([sehf; sd18o],1);
aves = aves./nanstd(aves);
figure; plot(cs,aves)
xlabel('Age (Ma)'); ylabel('Average Slope'); set(gca,'xdir','reverse');
formatfigure;

% Slope of average
saverage = saverage./nanstd(saverage);
figure; plot(cs,saverage)
xlabel('Age (Ma)'); ylabel('Slope of Average'); set(gca,'xdir','reverse');
formatfigure;

% Covariance
figure; plot(cs,covarianceHfO)
xlabel('Age (Ma)'); ylabel('Covariance'); set(gca,'xdir','reverse');
covslopeneg = covarianceHfO;
covslopeneg(aves>0) = NaN;
hold on; plot(cs,covslopeneg,'r')
formatfigure;

% % Covariance(>0) * Slope
% slopecov = covarianceHfO.*saverage;
% slopecov(covarianceHfO<0)=0;
% figure; plot(cs,slopecov./nanstd(slopecov))
% xlabel('Age (Ma)'); ylabel('Slope * Covariance')
% set(gca,'xdir','reverse')
% formatfigure

% abs(Covariance) * Slope
slopecovabs = abs(covarianceHfO).*saverage;
figure; plot(cs,slopecovabs)
xlabel('Age (Ma)'); ylabel('Slope * abs(Covariance)'); set(gca,'xdir','reverse');
formatfigure;
slopecovneg = slopecovabs;
slopecovneg(covarianceHfO>0)=NaN;
hold on; plot(cs,slopecovneg,'r')


%% Add glaciations

hold on;
plot([717,717],get(gca,'ylim')) %Sturtian
plot([580,580],get(gca,'ylim')) %Gaskiers
plot([635,635],get(gca,'ylim')) %Marinoan

% yl = get(gca,'ylim');
% fill([340 340 260 260],[yl fliplr(yl)],'b','FaceAlpha',0.5) %% LPIA

% plot([2180,2180],get(gca,'ylim')) %Makganyene

% plot([2220,2220],get(gca,'ylim')) %Ongeluk formation, diamictite end, Kopp et al, 2005; Rasmussen et al, 2013
% plot([2308,2308],get(gca,'ylim')) %Gordon lake Makganyene begin, Kopp et al, 2005; Rasmussen et al, 2013
% plot([2415,2415],get(gca,'ylim')) %Koegas, Kirschvink et al, 1999


plot([2220,2220],get(gca,'ylim')) %Rasmussen et al., 2013
plot([2320,2320],get(gca,'ylim')) %Rasmussen et al., 2013
plot([2370,2370],get(gca,'ylim')) %Rasmussen et al., 2013

plot([2415,2415],get(gca,'ylim')) %Koegas, Kirschvink et al, 1999

plot([2940,2940],get(gca,'ylim')) %Pongola, von Brunn and Gold, 1993


%% Plot by continent


Elem='eHf_initial';
agemin=0;
agemax=4350;
nbins = 145;
binoverlap = 3;

figure; hold on;

test=~isnan(mczircon.(Elem)) & mczircon.Continent==2; % Eurasia
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
errorbar(c,m,2*e,'.')

test=~isnan(mczircon.(Elem)) & (mczircon.Continent==3 | mczircon.Continent==4); % North or South America
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
errorbar(c,m,2*e,'.')

test=~isnan(mczircon.(Elem)) & mczircon.Continent==5; % Australia
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
errorbar(c,m,2*e,'.')

legend('Eurasia','Americas','Australia')


test=~isnan(mczircon.(Elem)) & (mczircon.Continent==1 | mczircon.Continent==7); % Africa and Antarctica
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
errorbar(c,m,2*e,'.')

xlabel('Age (Ma)'); ylabel(Elem); xlim([agemin, agemax])
set(gca,'xdir','reverse');
xlim([0 4000]);
% ylim([-10 10]);
set(gca,'ytick',[-10 -5 0 5 10])
formatfigure;
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig,'ZirconEHfSignalByContinent.pdf')


%%

mean(mcigncn1.K2O(~isnan(mcigncn1.K2O)))

mean(mcigncn1.Age(~isnan(mcigncn1.Age)))
