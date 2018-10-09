% mctest.m
% Produce the the (~10^6 row) monte carlo table for quick examination
% with plotmcvariable(s).m
% Uses per-sample relative errors
%%

% Elements to include when resampling
simitems={'Age';'d18O';'Hf176_Hf177';'Lu176_Hf177';'eHf_initial';'Discordance';'Continent';};

% Check that age data is present
agecol=strcmp(simitems,'Age');
if ~any(agecol); error('Error: missing age variable'); end;
if sum(agecol)>1; error('Error: multiple age variables'); end;

% Create data and uncertainty matrices
data=zeros(length(zircon.Age),length(simitems));
uncertainty=zeros(size(data));
for i=1:length(simitems)
    % Data matrix
    data(:,i)=zircon.(simitems{i});
    
    % Uncertainty matrix. Note that this works regardless of whether the
    % relative uncertainty is a vector or a scalar
    uncertainty(:,i)=zircon.err2srel.(simitems{i});
    
    % For any NaN uncertainties that aren't NaN data, set relative uncertainty to global average
    test=isnan(uncertainty(:,i))&~isnan(data(:,i));
    uncertainty(test,i)=mean(uncertainty(~isnan(uncertainty(:,i)),i));
end

% Set minimum absolute age uncertainty
MinAbsAgeUncert=50; 
test=uncertainty(:,agecol).*data(:,agecol) < MinAbsAgeUncert | isnan(uncertainty(:,agecol));
uncertainty(test,agecol)=MinAbsAgeUncert./data(test,agecol);


% Screen out any unwanted data
test=~isnan(zircon.Age); % No point resampling by age if there isn't age data
data=data(test,:);
uncertainty=uncertainty(test,:);

%% Produce sample weights

tic;
fprintf('Calculating sample weights: ')
% if isfield(zircon,'k')
%     k=zircon.k;
% else
    k=invweightAge(zircon.Age(test));
    zircon.k=k;
% end

prob=1./((k.*median(5./k))+1);
toc


%% Run the monte carlo

% Number of rows to simulate
samplerows=1E7;

tic;
fprintf('Allocating output matrix: ')
% Generate output data matrix
mczircon.data=NaN(samplerows,size(data,2));
toc

tic;
fprintf('Resampling: ')
% Fill the output data matrix with resampled data
i=1;
while i<samplerows
    % select weighted sample of data
    ru=rand(length(prob),1);
    sdata=data(prob>ru,:);
    suncertainty=uncertainty(prob>ru,:);
    
    % Randomize data over uncertainty interval
    rn=randn(size(sdata));
    sdata=sdata+rn.*suncertainty./2.*sdata;
    
    if i+size(sdata,1)-1<=samplerows
        mczircon.data(i:i+size(sdata,1)-1,:)=sdata;
    else
        mczircon.data(i:end,:)=sdata(1:samplerows-i+1,:);
    end
    i=i+size(sdata,1);
    
end
toc

% Add elements and separate data matrix into separate variables within the struct
mczircon.elements=simitems;
mczircon=elementify(mczircon);

% Recalculate eHf based on resampled Lu and Hf ratios
mczircon.eHf_initial=eHf(mczircon.Hf176_Hf177,mczircon.Lu176_Hf177,mczircon.Age);

% Save results
save mczircon mczircon

