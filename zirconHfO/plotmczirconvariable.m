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

test=~isnan(mczircon.(Elem));
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
figure; errorbar(c,m,2*e,'.')
xlabel('Age (Ma)'); ylabel(Elem); xlim([agemin, agemax])
set(gca,'xdir','reverse');
xlim([0 4500]);
ylim([-10 10]);
set(gca,'ytick',[-10 -5 0 5 10])
formatfigure

%% Plot Lu176/Hf177

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
formatfigure;

%% Plot d18O

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
formatfigure;


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

