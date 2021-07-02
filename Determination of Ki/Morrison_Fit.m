% Global Morrison Fit for Determination of Ki  with Statistical Packages
%
% Created: January 17, 2019
% Created by: Abby Harthorn, vandu054@umn.edu
% Corresponding author: Benjamin Hackel, hackel@umn.edu

% Description:
% Raw absorbance data from inhibition assay was saved to an Excel file and
% imported to MATLAB via xlsread. Excel files are formatted in six columns
% of background (no enzyme) at various inhibitor concentrations, followed
% by six columns of CA activity in presence of various inhibitor
% concentrations. Therefore, 12 data columns per trial, with at least 3
% replicates. Absorbance was read every minute for 30 minutes.
%
% Conditions: 25 nM CA-II; 100 nM CA-IX; 
% [I] = 200 nM, 100 nM, 50 nM, 10 nM, 1nM, 0nM

% Data file: FnII-28-7_Raw_Data includes raw data for AAZ and FnII.28.7-AAZ
% Raw data for other PriSMs were not included due to change in naming
% scheme. Additional data available upon request.

clear;
clc;
close all;

% Data file: FnII-28-7_Raw_Data
% sheet1 = AAZ CA-II
% sheet2 = AAZ CA-IX
% sheet3 = FnII.28.7-AAZ CA-II
% sheet4 = FnII.28.7-AAZ CA-IX

% To calculate all parameters for all samples at once, set A
% a corresponds to the Excel sheet (CA-II odd #s; CA-IX even #s)
for A=0:1
a = 2*A+1; 

% Total enzyme concentration (CA-II Et=25e-9; CA-IX Et=100e-9)
% Be sure to change to appropriate concentration 
Et = 25e-9; 

clear data vbgd vtarg

data = xlsread('FnII-28-7_Abs_Data',a); 

dpts = size(data,1);
trials = size(data,2)/12;

% If statement determines the number of replicates per sample
if trials == 2
    bgd = cat(3,data(:,1:6), data(:,13:18));
    targ = cat(3,data(:,7:12), data(:,19:24)); 
elseif trials == 3
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36));  
elseif trials == 4
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30), data(:,37:42));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36), data(:,43:48));
elseif trials == 5
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30), data(:,37:42), data(:,49:54));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36), data(:,43:48), data(:,55:60));        
elseif trials == 6 
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30), data(:,37:42), data(:,49:54), data(:,61:66));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36), data(:,43:48), data(:,55:60), data(:,67:72));
elseif trials == 7
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30), data(:,37:42), data(:,49:54), data(:,61:66),data(:,73:78));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36), data(:,43:48), data(:,55:60), data(:,67:72), data(:,79:84));
elseif trials == 8
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30), data(:,37:42), data(:,49:54), data(:,61:66),data(:,73:78),data(:,85:90));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36), data(:,43:48), data(:,55:60), data(:,67:72), data(:,79:84),data(:,91:96));
elseif trials == 9
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30), data(:,37:42), data(:,49:54), data(:,61:66),data(:,73:78),data(:,85:90),data(:,97:102));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36), data(:,43:48), data(:,55:60), data(:,67:72), data(:,79:84),data(:,91:96),data(:,103:108));
elseif trials == 10
    bgd = cat(3,data(:,1:6), data(:,13:18), data(:,25:30), data(:,37:42), data(:,49:54), data(:,61:66),data(:,73:78),data(:,85:90),data(:,97:102),data(:,109:114));
    targ = cat(3,data(:,7:12), data(:,19:24), data(:,31:36), data(:,43:48), data(:,55:60), data(:,67:72), data(:,79:84),data(:,91:96),data(:,103:108),data(:,115:120));
elseif trials > 10
    disp('MORE THAN 10 EXPERIMENTS!');
end
    

% Determine slope for each inhibitor condition at each sample
% Slope determination started at datapoint 3 because sometimes absorbance
% was noisy in the first couple datapoints.
dp1 = 3;

for j = 1:trials     %   each replicate (n = 3)
    for i=1:6   %   each concentration (6)
        xfit = dp1:dpts;
        ybgd = bgd(dp1:dpts,i,j)';
        ytarg = targ(dp1:dpts,i,j)';
        
        fit1= polyfit(xfit,ybgd,1);     % fits background to linear line y = mx + b -> displays fit1 = m b 
        vbgd(i,j)=fit1(1);              % this grabs the slope (1st polynomial coefficient), the velocity
        fit2 = polyfit(xfit,ytarg,1); % fits target to linear line
        vtarg(i,j) = fit2(1);         % grabs slope

        yfitbgd = polyval(fit1,dp1:dpts);       % data points for straight line, poly fit
        yresidbgd = ybgd-yfitbgd;               % residuals of fit (difference of data - fit)
        SSresid = sum(yresidbgd.^2);           % Sum of Squares of residuals
        SStotal = (length(ybgd)-1)*var(ybgd);      % SStot = sum of (ybgd - mean(ybgd))^2 -> E(X-u)^2 / N-1  = sigma^2 = variance
        rsqbgd(i,j) = 1-SSresid/SStotal; 
        
        yfittarg = polyval(fit2,dp1:dpts);
        yresidtarg = ytarg-yfittarg;
        SSresid2 = sum(yresidtarg.^2);
        SStotal2 = (length(ytarg)-1)*var(ytarg);
        rsqtarg(i,j) = 1-SSresid2/SStotal2;% 1 - (SSres/SStot) = R^2

    end
    
end

% Warning for bad slope fit.
if min(min([rsqbgd rsqtarg])) < .95
    disp('bad slope fit!!!!!')
end

%% Normalize to background

% Background slopes were averaged across all inhibitor concentrations, as
% inhibitor did not seem to affect 4-NPA hydrolysis
vbgdm = mean(vbgd(1:6,:)); 
stdbgd = std(vbgd(1:6,:)); % standard deviation

% Warning for high standard deviation of background slope values
if max(stdbgd) > 0.0006
    disp('check background values, high standard deviation');
end

% CA activity adjusted for background signal
vtarg_adj = vtarg - repmat(vbgdm,6,1); 

% Because Carbonic Anhydrase had different maximum activity on different
% days, values were normalized using the slope from enzyme only (no
% inhibitor)
vtnorm = vtarg_adj./vtarg_adj(6,:); 

% Reorganized matrix to single column vector
if size(vtnorm,2) == 2
    vtnorm2 = [vtnorm(:,1);vtnorm(:,2)];
elseif size(vtnorm,2) == 3
        vtnorm2 = [vtnorm(:,1);vtnorm(:,2);vtnorm(:,3)];
elseif size(vtnorm,2) == 4
        vtnorm2 = [vtnorm(:,1);vtnorm(:,2);vtnorm(:,3);vtnorm(:,4)];
elseif size(vtnorm,2) == 5
    vtnorm2 = [vtnorm(:,1);vtnorm(:,2);vtnorm(:,3);vtnorm(:,4);vtnorm(:,5)];
elseif size(vtnorm,2) == 6
    vtnorm2 = [vtnorm(:,1);vtnorm(:,2);vtnorm(:,3);vtnorm(:,4);vtnorm(:,5);vtnorm(:,6)];
elseif size(vtnorm,2) == 7
    vtnorm2 = [vtnorm(:,1);vtnorm(:,2);vtnorm(:,3);vtnorm(:,4);vtnorm(:,5);vtnorm(:,6);vtnorm(:,7)];
end

% Inhibitor concentrations used
It = [200e-9, 100e-9, 50e-9, 10e-9, 1e-9, 0e-9]';

% Morrison Equation Model for Tight Binding
% Beta0 is inital guess for Ki and Maximum Velocity/Activity
model = @(b,It) b(2)*(1-(((Et + It + b(1)) - (sqrt((Et + It + b(1)).^2 - 4*Et*It)))/(2 * Et)));
beta0 = [1e-9,1];

X = repmat(It,size(vtnorm,2),1);
Y = vtnorm2;

% Optimizatin options for nlinfit
options = statset('TolFun',1e-15,'TolX',1e-15,'MaxIter',1e4','Display','iter');

% Nonlinear Regression
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(X,Y,model,beta0,options);

% Determination of Confidence Interval and Standard Error
betaCovariance = inv(J.'*J)*var(R);
[ci,se] = nlparci2(beta,R,'jacobian',J,'alpha',0.32); %0.32

% Option to display values for Ki and maximum velocity
disp('beta'); disp(a); disp(beta)
disp('confidence interval'); disp(ci)


%%

% Values used for graphing model values
x = logspace(-10, -5,100);

% Values used for graphing inhibitor concentration values on log scale
It2 = [200e-9, 100e-9, 50e-9, 10e-9, 1e-9, 0.00001e-9]';

figure;           
semilogx(x,model(beta,x)/max(model(beta,x)),'Linewidth',2.0)
hold all
scatter(repmat(It,size(vtnorm,2),1),vtnorm2/max(model(beta,x)))
xlabel('[I] mol/L')
ylabel('reaction velocity')
title({''});
legend('fit','normalized data')
set(gca,'linewidth',2.0,'FontSize',18.0)
set(legend('fit','normalized data'),'Location','southwest');
ylim([-.1,1.1])
xlim([1e-10 5e-7])

% Values for all inhibitors can be stored in following matrices. Variables
% started with 'a_' to be easily seen at the top of the Workspace.

% Apparent Ki
a_k(A+1,:) = beta(1);
% Mean velocity at each [inhibitor]
a_vel(A+1,:) = mean(vtnorm,2)./max(model(beta,x)); %#ok<*SAGROW>
% Model fit values
a_fit(A+1,:) = model(beta,x)./max(model(beta,x));
% Standard error at each velocity
a_SE(A+1,:) = std(vtnorm'/max(model(beta,x)),1)/sqrt(trials);
% Confidence interval of Ki
a_CI(A+1,:) = ci(1,:);
% Standard deviation at each velocity
a_SD(A+1,:) = std(vtnorm'/max(model(beta,x)),1);
% Ki value; Maximum velocity value
a_param(A+1,:) = beta;
% Velocity of each slope, non-averaged
a_vpts(A+1,1:length(vtnorm2')) = vtnorm2';
% Residuals
a_R(A+1,:) = sum(R.^2);


end

