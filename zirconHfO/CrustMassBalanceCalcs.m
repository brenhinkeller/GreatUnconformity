if ~exist('mcigncn1','var')
    load mcigncn1
end
if ~exist('igncn1','var')
    load igncn1
end

if ~exist('mczircon','var')
    load mczircon
end
if ~exist('zircon','var')
    load zircon
end



% Recalculate initial zircon Hf composition;
[mczircon.eHf_initial, mczircon.Hf176_Hf177_initial]=eHf(mczircon.Hf176_Hf177,mczircon.Lu176_Hf177,mczircon.Age);


% Estimate initial magmatic Hf isotope composition;
[c,m,e] = bin(mczircon.Age,mczircon.Hf176_Hf177_initial,0,4500,1,45);
figure; errorbar(c,m,2*e); xlabel('Age (Ma)'); ylabel('Hf176/Hf177 initial');
mcigncn1.Hf176_Hf177_initial_approx = interp1(c,m,mcigncn1.Age);


%%

lambda = 1.867*10^-11; % Lutetium decay constant (Soderlund et al., 2004)
CHUR_Hf176_Hf177 = 0.282785; % Present-day CHUR Hf ratio (Bouvier et al., 2008)
CHUR_Lu176_Hf177 = 0.0336; %Present-day CHUR Lu/Hf (Bouvier et al., 2008)

% Calculate average crustal Hf isotope composition at all subsequent times
crust.date = (4000:-1:0)';

crust.Hf176_Hf177_t = NaN(size(crust.date));
crust.eHf_t = NaN(size(crust.date));

%%%%%%%%%%%%%%%%%%
mcigncn1.Lu176_Hf177 = mcigncn1.Lu*0.02599 ./ mcigncn1.Hf*0.186; 
mcigncn1.Lu176_Hf177(mcigncn1.Hf==0)=NaN;
%%%%%%%%%%%%%%%%%%

for i=1:length(crust.date)
    
    dt = mcigncn1.Age-crust.date(i);
    
    Hf176_Hf177_t = mcigncn1.Hf176_Hf177_initial_approx + mcigncn1.Lu176_Hf177.*(exp(mcigncn1.Age *10^6*lambda) - exp(crust.date(i) *10^6*lambda));
    
    crust.Hf176_Hf177_t(i) = nanmean(Hf176_Hf177_t(dt>0));
    
    % Calculate CHUR Hf ratio at time t
    CHUR_Hf176_Hf177_t = CHUR_Hf176_Hf177 - CHUR_Lu176_Hf177.*(exp(crust.date(i)*10^6*lambda) - 1);

    % Calculate corresponding epsilon Hf
    crust.eHf_t(i)=(crust.Hf176_Hf177_t(i) ./ CHUR_Hf176_Hf177_t - 1) .* 10^4;
    
end


figure; plot(crust.date, crust.eHf_t)
hold on; plot(crust.date,(4500-crust.date)*17/4500)
xlabel('Age (Ma)'); ylabel('Average crustal eHf');
set(gca,'xdir','reverse')



%%

mczircon.ContemporaryAverageCrust_eHf = interp1(crust.date,crust.eHf_t,mczircon.Age);
mczircon.ContemporaryDepletedMantle_eHf = (4500-mczircon.Age)*17/4500;
mczircon.HfCrustIncorporationPercentage = (mczircon.ContemporaryDepletedMantle_eHf-mczircon.eHf_initial)./(mczircon.ContemporaryDepletedMantle_eHf-mczircon.ContemporaryAverageCrust_eHf).*100;

Elem='HfCrustIncorporationPercentage';
agemin=0;
agemax=4350;
nbins = 145;
binoverlap = 3;

test=~isnan(mczircon.(Elem));
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
figure; yyaxis left
errorbar(c,m,2*e,'.')
xlabel('Age (Ma)'); ylabel('Percent average crust incorporated in new magmas');
set(gca,'xdir','reverse');
xlim([0 4500]);
ylim([0 100]);

yyaxis right
ylim([0 100]);
set(gca,'ydir','reverse');
ylabel('Percent mantle input')

fig = gcf;
fig.PaperSize = [fig.PaperPosition(3) fig.PaperPosition(4)];

saveas(fig,'PercentCrustMantleInNewMagmas.pdf')

