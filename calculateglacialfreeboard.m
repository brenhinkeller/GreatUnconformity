function [Freeboard] = calculateglacialfreeboard(age,dt,geotherm,yLithosphere,iceThickness,tGlaciation,totalExhumation,exhumationFraction,sedAccumulation,isostaticFraction,rho_mantle,rho_lithosphere,rho_water)
    K_m = 1/1000^2*31556952; % m^2/yr. Diffusivity (1 mm^2/s)
    alpha_m = 3*10^-5; % 1/K. Coeff. of thermal expansion

    % Sediment accumulation and consequent freeboard
        sedFreeboard = sedAccumulation * 1000/isostaticFraction; % m

    % Glacial exhumation and consequent freeboard
        glacialExhumation = zeros(size(age));
        iceFreeboard = zeros(size(age));
        for i=1:size(tGlaciation,1)
            glacialExhumation(age<tGlaciation(i,2)) = glacialExhumation(age<tGlaciation(i,2)) + totalExhumation*exhumationFraction(i); % km
            iceFreeboard(age<tGlaciation(i,1) & age>=tGlaciation(i,2)) = NaN; % This gets complicated -70 + iceThickness*1000*0.3/0.7; % m
        end
        glacialFreeboard = -glacialExhumation*1000/isostaticFraction; % m

    % Immediate isostatic thermal difference
        Delta_t = geotherm*glacialExhumation;
        PThermalIsostatic = yLithosphere*1000*rho_lithosphere*(Delta_t/2*alpha_m);
                
    % Gradual lithospheric cooling and subsidence
        
        % Long-term lithospheric subsidenc
        % PResidualSubsidence = 2*rho_lithosphere*alpha_m*(1300-0) * (sqrt(K_m*(3000+t)*1e6/pi) - sqrt(K_m*3000*1e6/pi));
        PResidualSubsidence = zeros(size(age));
        m = 20;
        tm = repmat((0:dt:2000)',1,m+1);
        mm = repmat(0:m,size(tm,1),1);
        CratonicPThermalSubsidence = rho_lithosphere.*alpha_m.*(1300-0).*yLithosphere.*1000.*...
            (1/2 - 4/pi^2*sum(exp(-(K_m*(1+2*mm).^2*pi^2.*tm*1e6)/(yLithosphere*1000)^2)./(1+2*mm).^2,2));
        PResidualSubsidence = CratonicPThermalSubsidence(end-max(age)/dt:end)-CratonicPThermalSubsidence(end-max(age)/dt);
        
        % Subsidence in response to glacial exhumation
        % PThermalAdjustment = 2*rho_lithosphere*alpha_m*(Delta_t-0) * sqrt(K_m*t*1e6/pi);
        PThermalSubsidence = zeros(size(age));
        for i=1:size(tGlaciation,1)
            m = 20;
            tm = repmat((0:dt:tGlaciation(i,2))',1,m+1);
            mm = repmat(0:m,size(tm,1),1);
            PThermalSubsidence(age<=tGlaciation(i,2)) = PThermalSubsidence(age<=tGlaciation(i,2)) + ...
                rho_lithosphere.*alpha_m.*(geotherm*totalExhumation*exhumationFraction(i)-0).*yLithosphere.*1000.*...
                (1/2 - 4/pi^2*sum(exp(-(K_m*(1+2*mm).^2*pi^2.*tm*1e6)/(yLithosphere*1000)^2)./(1+2*mm).^2,2));
        end

        thermalFreeboard = (PThermalIsostatic-PResidualSubsidence-PThermalSubsidence)/(rho_mantle-rho_water);

        Freeboard = iceFreeboard + glacialFreeboard+thermalFreeboard+sedFreeboard;

end
