% %% Assign continents, Greenland as separate continent
% 
% load zircon
% load zircontext
% 
% empties = cellfun(@isempty,zircontext.Continent);
% 
% zircontext.Continent(empties) = cellfun(@(x) '',zircontext.Continent(empties),'UniformOutput',false);
% 
% zircon.Continent=NaN(size(zircontext.Continent));
% 
% % Africa
% parts = {'Africa','Afrique','Barberton','Swaziland','Zimbabwe'};
% for i=1:length(parts)
%     test = contains(zircontext.Continent,parts{i});
% %     test = strcmp(parts{i},zircontext.Continent);
%     zircon.Continent(test) = 1;
% end
% 
% % Eurasia
% parts = {'Asia','Europe','China','Chine','Eurasia','Norway','Finland','Russi','Siberi','Inde','India'};
% for i=1:length(parts)
%     test = contains(zircontext.Continent,parts{i});
% %     test = strcmp(parts{i},zircontext.Continent);
%     zircon.Continent(test) = 2;
% end
% 
% % North America
% parts = {'North America','N.Am','Quebec','USA','Mississippi','Canada','Acasta'};
% for i=1:length(parts)
%     test = contains(zircontext.Continent,parts{i});
% %     test = strcmp(parts{i},zircontext.Continent);
%     zircon.Continent(test) = 3;
% end
% 
% % South America
% parts = {'South America','S.Am','Brazil'};
% for i=1:length(parts)
%     test = contains(zircontext.Continent,parts{i});
% %     test = strcmp(parts{i},zircontext.Continent);
%     zircon.Continent(test) = 4;
% end
% 
% % Australia
% parts = {'Australi','Jack Hills','Pilbara','Yilgarn'};
% for i=1:length(parts)
%     test = contains(zircontext.Continent,parts{i});
% %     test = strcmp(parts{i},zircontext.Continent);
%     zircon.Continent(test) = 5;
% end
% 
% % Greenland
% parts = {'Greenland','Greenland','groenland'};
% for i=1:length(parts)
%     test = contains(zircontext.Continent,parts{i});
% %     test = strcmp(parts{i},zircontext.Continent);
%     zircon.Continent(test) = 6;
% end
% 
% % Antarctica
% parts = {'Antarcti'};
% for i=1:length(parts)
%     test = contains(zircontext.Continent,parts{i});
% %     test = strcmp(parts{i},zircontext.Continent);
%     zircon.Continent(test) = 7;
% end
% 
% dt=200;
% t=(0:dt:4400)';
% 
% Africa=histcounts(zircon.Age(zircon.Continent == 1),t)';
% Eurasia=histcounts(zircon.Age(zircon.Continent == 2),t)';
% NorthAmerica=histcounts(zircon.Age(zircon.Continent == 3),t)';
% SouthAmerica=histcounts(zircon.Age(zircon.Continent == 4),t)';
% Australia=histcounts(zircon.Age(zircon.Continent == 5),t)';
% Greenland=histcounts(zircon.Age(zircon.Continent == 6),t)';
% Antarctica=histcounts(zircon.Age(zircon.Continent == 7),t)';
% 
% 
% total=Africa+Eurasia+NorthAmerica+SouthAmerica+Australia+Greenland+Antarctica;
% 
% 
% figure;
% bar(center(t),[Africa Eurasia NorthAmerica SouthAmerica Australia Greenland Antarctica]./repmat(total,1,7)*100,'BarLayout','stacked','BarWidth',1)
% 
% ylim([1,100])
% legend('Africa','Eurasia','North America','South America','Australia','Greenland','Antarctica')
% set(gca,'XTick',0:1000:4000)
% set(gca,'YTick',0:20:100)
% xlabel('Age (Ma)');
% ylabel('Percent of zircon record');
% set(gca,'Xdir','reverse');
% formatfigure;

%% Assign continents, Greenland as part of North America

load zircon
load zircontext

empties = cellfun(@isempty,zircontext.Continent);

zircontext.Continent(empties) = cellfun(@(x) '',zircontext.Continent(empties),'UniformOutput',false);
zircon.Continent=NaN(size(zircontext.Continent));
zircon.err2srel.Continent = 0; % No uncertainty for resampling

% Africa
parts = {'Africa','Afrique','Barberton','Swaziland','Zimbabwe'};
for i=1:length(parts)
    test = contains(zircontext.Continent,parts{i});
    zircon.Continent(test) = 1;
end

% Eurasia
parts = {'Asia','Europe','China','Chine','Eurasia','Norway','Finland','Russi','Siberi','Inde','India'};
for i=1:length(parts)
    test = contains(zircontext.Continent,parts{i});
    zircon.Continent(test) = 2;
end

% North America
parts = {'North America','N.Am','Quebec','USA','Mississippi','Canada','Acasta','Greenland','Greenland','groenland'};
for i=1:length(parts)
    test = contains(zircontext.Continent,parts{i});
    zircon.Continent(test) = 3;
end

% South America
parts = {'South America','S.Am','Brazil'};
for i=1:length(parts)
    test = contains(zircontext.Continent,parts{i});
    zircon.Continent(test) = 4;
end

% Australia
parts = {'Australi','Jack Hills','Pilbara','Yilgarn'};
for i=1:length(parts)
    test = contains(zircontext.Continent,parts{i});
    zircon.Continent(test) = 5;
end

% Antarctica
parts = {'Antarcti'};
for i=1:length(parts)
    test = contains(zircontext.Continent,parts{i});
    zircon.Continent(test) = 7;
end



%% Plot raw (non-MC'd)abundance
dt=30; % Age step to plot with 
t=(0:dt:4350)';

Africa=histcounts(zircon.Age(zircon.Continent == 1),t)';
Eurasia=histcounts(zircon.Age(zircon.Continent == 2),t)';
NorthAmerica=histcounts(zircon.Age(zircon.Continent == 3),t)';
SouthAmerica=histcounts(zircon.Age(zircon.Continent == 4),t)';
Australia=histcounts(zircon.Age(zircon.Continent == 5),t)';
Antarctica=histcounts(zircon.Age(zircon.Continent == 7),t)';

total=Africa+Eurasia+NorthAmerica+SouthAmerica+Australia+Antarctica;


figure;
bar(center(t),[Africa Eurasia NorthAmerica SouthAmerica Australia Antarctica]./repmat(total,1,6)*100,'BarLayout','stacked','BarWidth',1)

ylim([1,100])
xlim([0,4500])
legend('Africa','Eurasia','North America','South America','Australia','Antarctica')
set(gca,'XTick',0:1000:4000)
set(gca,'YTick',0:20:100)
xlabel('Age (Ma)');
ylabel('Percent of zircon record');
set(gca,'Xdir','reverse');
formatfigure;


%% Plot Monte-Carlo'd abundances
load mczircon
dt=30; % Age step to plot with 
t=(0:dt:4350)';

Africa=histcounts(mczircon.Age(mczircon.Continent == 1),t)';
Eurasia=histcounts(mczircon.Age(mczircon.Continent == 2),t)';
NorthAmerica=histcounts(mczircon.Age(mczircon.Continent == 3),t)';
SouthAmerica=histcounts(mczircon.Age(mczircon.Continent == 4),t)';
Australia=histcounts(mczircon.Age(mczircon.Continent == 5),t)';
Antarctica=histcounts(mczircon.Age(mczircon.Continent == 7),t)';

total=Africa+Eurasia+NorthAmerica+SouthAmerica+Australia+Antarctica;


figure;
bar(center(t),[Africa Eurasia NorthAmerica SouthAmerica Australia Antarctica]./repmat(total,1,6)*100,'BarLayout','stacked','BarWidth',1)

ylim([1,100])
xlim([0,4500])
legend('Africa','Eurasia','North America','South America','Australia','Antarctica')
set(gca,'XTick',0:1000:4000)
set(gca,'YTick',0:20:100)
xlabel('Age (Ma)');
ylabel('Percent of zircon record');
set(gca,'Xdir','reverse');
formatfigure;
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig,['ZirconContinentDistribution' num2str(dt) '.pdf'])

