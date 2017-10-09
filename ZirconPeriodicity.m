load zircon
load mczircon
load voice

agemin=0;
agemax=4350;
% nbins = 145;
% binoverlap = 3;
nbins = 290;
binoverlap = 0;

% Calculate zircon eHf timeseries
Elem='eHf_initial';
test=~isnan(mczircon.(Elem));
[c,ehf,ehf_sigma]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);

% Calculate ziron eHf timeseries
Elem='d18O';
test=~isnan(mczircon.(Elem));
[c,d18o,d18o_sigma]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);

data.BinCenterAge = c;
data.eHf = ehf;
data.ehf_sigma = ehf_sigma;
data.d18O = d18o;
data.d18O_sigma = d18o_sigma;
data.elements = fieldnames(data);
exportdataset(data,'ZirconTimeseriesData.csv',',');

% De-trend and unit variance
d18o = detrend(d18o);
d18o = d18o/std(d18o);
ehf = detrend(ehf);
ehf = -ehf/std(ehf); % Negative since lower eHf means more crust remelting
figure; plot(c,ehf)
figure; plot(c,d18o)



%%
freq = 0.25:0.01:20; % Cycles/Gyr

[pxx1,f] = periodogram(ehf,[],freq,1000/15);
figure; hold on; plot(1000./f,pxx1)

[pxx2,f] = periodogram(d18o,[],freq,1000/15);
plot(1000./f,pxx2)
legend('eHf','d18O')
xlabel('Period (Myr)'); ylabel('Spectral power');
xlim([0 1009])
formatfigure;

%%

freq = 0.25:0.01:20; % Cycles/Gyr

[pxx1,f] = periodogram(ehf,[],freq,1000/15);
figure(1); hold on; plot(1000./f,pxx1)
% figure(2); hold on; plot(f,pxx1)

[pxx2,f] = periodogram(d18o,[],freq,1000/15);
figure(1); plot(1000./f,pxx2)
legend('eHf','d18O')
% figure(2); plot(f,pxx2)

% spacing = 30;
% edges = 100:spacing:4100;
% Age = center(edges);
% N = histcounts(zircon.Age,edges);
spacing = 1;
[N,Age] = ksdensity(zircon.Age,0:1:4500,'bandwidth',10);
figure(3); hold on; plot(Age,N./max(N))
N = detrend(N);
N = N./nanstd(N);
[pxx3,f] = periodogram(N,[],freq,1000/spacing);
figure(1); plot(1000./f,pxx3)
% figure(2); plot(f,pxx3)


[N,Age] = ksdensity(voice.Best_Age,0:1:4500,'bandwidth',10);
figure(3); hold on; plot(Age,N./max(N))
N = detrend(N);
N = N./nanstd(N);
[pxx4,f] = periodogram(N,[],freq,1000/1);
figure(1); plot(1000./f,pxx4)
xlabel('Period (Myr)'); ylabel('Spectral power');
legend('eHf','d18O','Zircon abundance','Voice zircon abundance')
xlim([0 1800])
% figure(2); plot(f,pxx4)
% xlabel('Frequency (Gyr^-^1)'); ylabel('Spectral power');
% legend('eHf','d18O','Zircon abundance','Voice zircon abundance')
% xlim([0 5])

%%
figure; plot(1000./f,nanmean([pxx1; pxx2; pxx3]))
xlabel('Period (Myr)'); ylabel('Average spectral power');
xlim([0 1800])
formatfigure

