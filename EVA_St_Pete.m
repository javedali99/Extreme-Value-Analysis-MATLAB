load WaterLevel_St_Pete_hourly.mat


%% Obtain annual maxima (AMAX) and MSL data from hourly record
% change datum format
Tmat = datevec(Time);

% get (unique) years for which we have data
years = unique(Tmat(:,1));

% get number of years
Nyear = length(years); 

% Use loop to identify annual maxima
for ii = 1:Nyear
    % find all values from year ii
    V = find(Tmat(:,1) == years(ii));
    AMAX(ii,1) = nanmax(WaterLevel(V));
    Tamax(ii,1) = years(ii);
    % MSL(ii,1) = nanmean(h_dat(V);  
end

clear years V ii h_dat Tmat T

% plot AMAX time series
figure
plot(Tamax,AMAX,'color',[0 0 0],'marker','o','markerfacecolor',[1 1 1],'linewidth',2)
xlabel('Time (yrs)')
ylabel('Level (mMSL)')

%% Remove the linear trend
Para = polyfit(Tamax,AMAX,1);
LinTrend = Para(1)*Tamax+Para(2);

% We can use rergress to check for significance!
[b,bint] = regress(AMAX,[ones(66,1) Tamax]);

% plot AMAX time series with linear trend
figure
plot(Tamax,AMAX,'color',[0 0 0],'marker','o','markerfacecolor',[1 1 1],'linewidth',2)
hold on
plot(Tamax,LinTrend,'color','r') 
xlabel('Time (yrs)')
ylabel('Level (mMSL)')

% get difference between last value in LinTrend and the other values
Diff = LinTrend(end,1) - LinTrend;
AMAXdt = AMAX+Diff;

% test if new trend is zero
Para2 = polyfit(Tamax,AMAXdt,1);
LinTrend2 = Para2(1)*Tamax+Para2(2);

% plot AMAX time series and detrended AMAX 
figure
plot(Tamax,AMAX,'color',[0 0 0],'marker','o','markerfacecolor',[1 1 1],'linewidth',2)
hold on
plot(Tamax,LinTrend,'color',[0 0 0])
plot(Tamax,AMAXdt,'color','r','marker','o','markerfacecolor',[1 1 1],'linewidth',2)
plot(Tamax,LinTrend2,'color','r')
xlabel('Time (yrs)')
ylabel('Level (mMSL)')

%% Detrend by subtracting mean sea level
%AMAXdt2 = AMAX-MSL;


%% Get Weibull plotting positions
% Check this: https://www.mathworks.com/help/stats/wblplot.html
WBL = [1:Nyear]'/(Nyear+1); % it gives non-exceedance probability
WBLrp = 1./(1-WBL); %return period
So = sort(AMAXdt); %sorting the values


figure
semilogx(WBLrp,So,'marker','+','color',[0 0 0],'linestyle','none')
xlabel('Return period (yrs)')
ylabel('Level (mMSL)')

%% Get GEV parameters and values for selected return periods
% fit GEV on detrended annual maxima series and get parameters
GEVpara = gevfit(AMAXdt); %default method: Maximum Likelihood Method

% define return periods we are interested in
rp = [1.01 1.1 1.25 1.5 2 5 10 50 100 250 500 1000];
Pu = 1-1./rp; %non-exceedance probability

% use inverse GEV to obtain return levels associated with the selected
% return periods
RL = gevinv(Pu,GEVpara(1),GEVpara(2),GEVpara(3)); % RL = return level

%% Alternatively, get return period for given level (1.9 m)
Non_exceedance = gevcdf(1.9, GEVpara(1),GEVpara(2), GEVpara(3));
Exceedance = 1-Non_exceedance;
Return_period = 1/Exceedance;

%% Get exceedance probabilities for 1.3m, 1.7m and 2.1m water levels
Non_exceedance1 = gevcdf(1.3, GEVpara(1),GEVpara(2), GEVpara(3));
Exceedance1 = 1-Non_exceedance1;
Return_period1 = 1/Exceedance1;

Non_exceedance2 = gevcdf(1.7, GEVpara(1),GEVpara(2), GEVpara(3));
Exceedance2 = 1-Non_exceedance2;
Return_period2 = 1/Exceedance2;

Non_exceedance3 = gevcdf(2.1, GEVpara(1),GEVpara(2), GEVpara(3));
Exceedance3 = 1-Non_exceedance3;
Return_period3 = 1/Exceedance3;

%% Return period/level plot
figure
semilogx(WBLrp,So,'marker','+','color',[0 0 0],'linestyle','none')
hold on
semilogx(rp,RL,'color','r','linewidth',2)
grid on
xlabel('Return period (yrs)')
ylabel('Level (mMSL)')


%% Assume sea level rise

% Add 20 cm to location paramater
GEVloc20 = GEVpara(3)+0.2;

RL20 = gevinv(Pu,GEVpara(1),GEVpara(2),GEVloc20);

% Plot
figure
semilogx(WBLrp,So,'marker','+','color',[0 0 0],'linestyle','none')
hold on
semilogx(rp,RL,'color','r','linewidth',2)
semilogx(rp,RL20,'color','b','linewidth',2)
grid on
xlabel('Return period (yrs)')
ylabel('Level (mMSL)')

