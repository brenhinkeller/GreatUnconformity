%% Estimate thermal and isostatic contributions

    K_m = 1/1000^2*31556952;	% m^2/yr. Diffusivity (1 mm^2/s)
    alpha_m = 3*10^-5;          % 1/K. Coeff. of thermal expansion
    areaLand = 148940000;       % km^2. Land area
    geotherm = 25;              % C/km

    rho_mantle = 3300;      % kg/m^3. Mantle density
    rho_lithosphere = 3250; % kg/m^3
    rho_crust = 2818;       % kg/m^3
    rho_water = 1000;       % kg/m^3. Water density

    yLithosphere = 133; % km
    yCrust = 33;        % km
    iceThickness = 2;   % km
    NormalizeSeaLevel = true;

    dt = 1;             % Myr
    t = (0:dt:800)';    % Myr
    age = flipud(t);    % Ma
    
    if ~exist('macrostrat','var'); load ../macrostrat/macrostrat.mat; end
    macrostratIndex = findclosest(flipud(t),macrostrat.Time-0.5); % Match macrostrat ages 

    sedRate = macrostrat.Volume(macrostratIndex)*6.07/1e6;  % Macrostrat sed rate
    sedRate = 0.9*ones(size(t));                            % Constant sed rate

    
    tGlaciation = [717.4, 660;     % Ma
                   640.7, 635.5;    % Ma
                   581, 580];	% Ma
               
    exhumationFraction = [1/3; 1/3; 1/3]; % Equal exhumation in each glaciation
%     exhumationFraction = tGlaciation(:,1)-tGlaciation(:,2); % Exhumation proportional to glacial duration
    exhumationFraction = exhumationFraction/sum(exhumationFraction); % Normalize

    
    
%% Vary magnitude of glacial exhumation
    exhumationRange = 1:2:9; % km
    yLithosphere = 150; % km
    
    sedAccumulation = cumsum(1e6*dt * sedRate) / areaLand; % km
    figure;
    isostaticFraction = (rho_crust-rho_water) / (rho_mantle - ((yLithosphere-yCrust)*rho_lithosphere + (yCrust)*rho_crust) / yLithosphere); % Crustal and lithospheric isostasy % Assume 1/2 of crust is generally covered by ocean
    for totalExhumation = exhumationRange;
        Freeboard = calculateglacialfreeboard(age,dt,geotherm,yLithosphere,iceThickness,tGlaciation,totalExhumation,exhumationFraction,sedAccumulation,isostaticFraction,rho_mantle,rho_lithosphere,rho_water);
        if NormalizeSeaLevel
            Freeboard = Freeboard-Freeboard(end); % Normalize to present-day sea level
        end
        hold on; plot(age,Freeboard)
    end
    legend(cellfun(@num2str, num2cell(exhumationRange),'UniformOutput',false));
    ylabel('Freeboard'); xlabel('Age (Ma)')
    set(gca,'XDir','Reverse')
%     yl = ylim;
%     ylim([-max(yl),max(yl)])
    formatfigure
    
%% Vary Lithospheric thickness
    yLithosphereRange = 50:100:350;
    totalExhumation = 3.35;
    totalExhumation = 4.5;

    sedAccumulation = cumsum(1e6*dt * sedRate) / areaLand; % km
    figure;
    for yLithosphere = yLithosphereRange;
        isostaticFraction = (rho_crust-rho_water) / (rho_mantle - ((yLithosphere-yCrust)*rho_lithosphere + (yCrust)*rho_crust) / yLithosphere); % Crustal and lithospheric isostasy % Assume 1/2 of crust is generally covered by ocean
        Freeboard = calculateglacialfreeboard(age,dt,geotherm,yLithosphere,iceThickness,tGlaciation,totalExhumation,exhumationFraction,sedAccumulation,isostaticFraction,rho_mantle,rho_lithosphere,rho_water);
        if NormalizeSeaLevel
            Freeboard = Freeboard-Freeboard(end); % Normalize to present-day sea level
        end
        hold on; plot(age,-Freeboard)
    end
    legend(cellfun(@num2str, num2cell(yLithosphereRange),'UniformOutput',false));
    ylabel('Sea Level'); xlabel('Age (Ma)')
    set(gca,'XDir','Reverse')
    yl = ylim;
    hold on; plot([541,541],ylim,'--k')
    ylim([-max(yl),max(yl)])
    
%% Pick a good pair
    yLithosphere = 133;
    
    sedRate = macrostrat.Volume(macrostratIndex)*6.07/1e6;  % Macrostrat sed rate
    totalExhumation = 3.35;
    
%     sedRate = 0.9*ones(size(t));                            % Constant sed rate
%     totalExhumation = 4.5;

    
    sedAccumulation = cumsum(1e6*dt * sedRate) / areaLand; % km
    isostaticFraction = (rho_crust-rho_water) / (rho_mantle - ((yLithosphere-yCrust)*rho_lithosphere + (yCrust)*rho_crust) / yLithosphere); % Crustal and lithospheric isostasy % Assume 1/2 of crust is generally covered by ocean
    Freeboard = calculateglacialfreeboard(age,dt,geotherm,yLithosphere,iceThickness,tGlaciation,totalExhumation,exhumationFraction,sedAccumulation,isostaticFraction,rho_mantle,rho_lithosphere,rho_water);
    if NormalizeSeaLevel
        Freeboard = Freeboard-Freeboard(end); % Normalize to present-day sea level
    end
    figure; hold on;
    plot(age,Freeboard)
    
    legend(num2str(yLithosphere));
    xlabel('Age (Ma)'); ylabel('Freeboard');
    set(gca,'XDir','Reverse')
    formatfigure;

    
%% Calculate fractional coverage relative to modern-day hypsometry

    if ~exist('etopoelev','var'); load etopoelev; end
    if ~exist('etoposmall','var'); etoposmall = imresize(etopoelev,1/8); end

    uncoveredFraction = NaN(size(t));
    for i=1:length(t)
        uncoveredFraction(i) = sum(etoposmall(:)>-Freeboard(i))/sum(etoposmall(:)>-40);
    end
    uncoveredFraction(isnan(Freeboard))=NaN;
    figure; hold on;
    plot(flipud(t),1-uncoveredFraction)
    set(gca,'Xdir','Reverse')
    xl = xlim;
    xlabel('Age (Ma)'); ylabel('Covered Fraction');
    plot(macrostrat.Time,macrostrat.Coverage)
    xlim(xl);
    legend('Model','Macrostrat')
    formatfigure;
    
    
%% Add Ronov coverage fraction
    
    RonovCoverage = importdataset('ronov_cont_flood.csv',',');
    
    t = contains(RonovCoverage.setting,'WHOLE CONTINENT');
    
    plot(RonovCoverage.mid_age(t),RonovCoverage.percent_flooded(t)/100,'.-')

    
    