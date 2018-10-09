load zircon

%% Plot O dataset

figure; plot(zircon.Age/1000, zircon.d18O,'.')
xlabel('Age (Ga)'); ylabel('\delta^{18}O')
set(gca,'xdir','reverse')
formatfigure;
ylim([0 14])

warning('off', 'MATLAB:print:FigureTooLargeForPage')
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig,'ZirconO.pdf')


%% Plot Hf dataset

figure; plot(zircon.Age/1000, zircon.eHf_initial,'.')
xlabel('Age (Ga)'); ylabel('\epsilonHf')
set(gca,'xdir','reverse')
formatfigure;
ylim([-60, 20])

% Add depleted mantle
hold on; plot(0:0.1:4.5,(4.5-(0:0.1:4.5))*17/4.5)
fig = gcf; fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];
saveas(fig,'ZirconHf.pdf')

