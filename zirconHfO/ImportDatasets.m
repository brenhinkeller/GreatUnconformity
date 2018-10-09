%% Load datasets
tic;

[TerraneChron, TerraneChronText]=importdataset('datasets/Belousova/ZirconHf.csv','|');

[UCLA, UCLAtext]=importdataset('datasets/UCLA/ZirconHf.csv','|');

[Dhuime, Dhuimetext]=importdataset('datasets/Dhuime/ZirconHf.csv','|');

[DhuimeHfO, DhuimeHfOtext]=importdataset('datasets/Dhuime/ZirconHfO.csv','|');

[PayneHfO, PayneHfOtext]=importdataset('datasets/Payne/ZirconHfO.csv','|');

[PayneO, PayneOtext]=importdataset('datasets/Payne/ZirconO.csv','|');

toc;

%% Combine datasets

zircon=concatenatedatasets(UCLA, TerraneChron);
zircon=concatenatedatasets(zircon, Dhuime);
zircon=concatenatedatasets(zircon, DhuimeHfO);
zircon=concatenatedatasets(zircon, PayneHfO);
zircon=concatenatedatasets(zircon, PayneO);
% save zircon zircon

zircontext=concatenatedatasets(UCLAtext, TerraneChronText);
zircontext=concatenatedatasets(zircontext, Dhuimetext);
zircontext=concatenatedatasets(zircontext, DhuimeHfOtext);
zircontext=concatenatedatasets(zircontext, PayneHfOtext);
zircontext=concatenatedatasets(zircontext, PayneOtext);
% save zircontext zircontext


%% Move errors to .err struct as relative uncertainties

for e=zircon.elements'
    % add default uncertainty of 1 per mil
    zircon.err2srel.(e{1})=0.001;
end


for elem = zircon.elements'
    e=elem{1};
    
    % Add _2sigma uncertainties
    tryname=[e '_2sigma'];
    if isfield(zircon,tryname)
        if length(zircon.err2srel.(e))==length(zircon.(tryname)) % Just add new data if variable already exists
            zircon.err2srel.(e)=nanmean([zircon.err2srel.(e), zircon.(tryname)./zircon.(e)],2);
        else
            zircon.err2srel.(e)=zircon.(tryname)./zircon.(e);
        end
        zircon=rmfield(zircon,tryname);
        zircon.err2srel=rmfield(zircon.err2srel,tryname);
        zircon.elements=zircon.elements(~strcmp(zircon.elements,tryname));
    end
    
    % Add _1sigma uncertainties
    tryname=[e '_1sigma'];
    if isfield(zircon,tryname)
        if length(zircon.err2srel.(e))==length(zircon.(tryname))
            zircon.err2srel.(e)=nanmean([zircon.err2srel.(e), 2.*zircon.(tryname)./zircon.(e)],2);
        else
            zircon.err2srel.(e)=2.*zircon.(tryname)./zircon.(e);
        end
        zircon=rmfield(zircon,tryname);
        zircon.err2srel=rmfield(zircon.err2srel,tryname);
        zircon.elements=zircon.elements(~strcmp(zircon.elements,tryname));
    end
end

for elem = zircon.elements'
    % Zero-valued uncertainties aren't meaningful; NaN out
    zircon.err2srel.(elem{1})(zircon.err2srel.(elem{1})==0)=NaN;
    zircon.err2srel.(elem{1})(zircon.err2srel.(elem{1})==Inf)=NaN;
    zircon.err2srel.(elem{1})(zircon.err2srel.(elem{1})==-Inf)=NaN;
    
    % Set unknown uncertainties based on the average uncertainty
    if size(zircon.err2srel.(elem{1})) == size(zircon.(elem{1}))
        test = isnan(zircon.err2srel.(elem{1})) & ~isnan(zircon.(elem{1}));
        zircon.err2srel.(elem{1})(test) = nanmedian(zircon.err2srel.(elem{1}));
    end
end


%% Recalculate epsilon Hf and eliminate major outliers;
zircon.eHf_initial=eHf(zircon.Hf176_Hf177,zircon.Lu176_Hf177,zircon.Age);

% test=zircon.eHf_initial > (24 - 20./4700.*zircon.Age);
% zircon.eHf_initial(test)=NaN;

save zircon zircon



%% Check how many rows overlap

% Just the essentials
checkelements={'Age';'Hf176_Hf177';'Lu176_Hf177';'d18O'};

% Round down to four sig figs to make sure different rounding lengths
% between databases aren't causing the samples to fail an equality test
testmatrix=zeros(length(zircon.Age),length(checkelements));
for i=1:length(checkelements)
    testmatrix(:,i)=roundsigfigs(zircon.(checkelements{i}),4);
end

% Set NaNs to zero for comparison purposes
testmatrix(isnan(testmatrix))=0;

% Test for uniqueness
uniques=0;
for i=1:length(testmatrix)
    if ~any(all(testmatrix(i+1:end,:)==repmat(testmatrix(i,:),length(testmatrix)-i,1),2))
        uniques=uniques+1;
    end
end

display(uniques);
