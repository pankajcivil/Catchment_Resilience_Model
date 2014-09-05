function monthlyLAI = sinLAI( LAImin, LAImax, LAIphase, t)
%sinMonthlyPrecip Summary of this function goes here
    
    % get within year time
    tmonth = t - floor(t);
    
    monthlyLAI(:,1) = LAImin + 0.5*(LAImax - LAImin).*(1-cos(2.*pi().*(tmonth+LAIphase)));
