%% Down payment analyzer
clc; clear; close all;

% This script determines how much value there is in putting a down payment
% on a home. Factors to consider include the home price, the down payment
% percentage, and how much value that down payment would accrue if instead
% we were to rent and let the down payment sit in the trust. For this
% analysis, I am using market values for interest rates (3.5APR), property 
% taxes in Nashville, average homeowners insurance, assumed 4% annual 
% housing market appreciation, and 7% annual trust appreciation. I am
% assuming that the trust appreciation accounts for investment taxes and
% financial managment fees. 

% Furthermore, I am assuming that the difference between monthly costs for
% owning (total mortgage payment, HO insurance, property taxes, etc.) and
% the monthly cost for renting (rent) is reinvested into the market
% (trust). I am also assuming that there is a 6% realtor fee that we pay
% up-front. From this analysis, it appears that renting is actually most
% cost-effective in the short-term (planning to sell within 10y). If
% buying, it is most cost effective to put as little down as possible. 
% For any down payment value, home ownership is cost effective over the 
% renting/investing strategy in the long term.

% This analysis is done for a house price of 750k, but the same takeaways
% hold for other house prices.


%% Script options

plot_mortgage = true;
plot_ownership_value = true;
plot_renting_value = true;
plot_difference = true;



%% Define variables

dps = [0.05, .1, .15, .2, .25];  % down payment percentages to consider
iter = 1;
for dp = dps
    a = 1.07;   % annual housing appreciation (to be compounded monthly)
    p = 650000; % house price
    d = p*dp;   % down payment amount
    f = 1.035;  % annual mortgage interest
    b = 1.08;   % annual trust appreciation
    r = 2300;   % monthly rent
    h = 0.03155 / 4;   % annual property taxes (3.155% on 25% of appraisal)
    N = 360;    % months in mortgage schedule
    z = 1200;    % misc monthly fees (homeowners insurance, lending fees, etc)
    

    %% Compute monthly mortgage payment

    m = p - d; % principle mortgage loan
    f_m = (f-1)/12; % monthly mortgage rate
    g = m * f_m*(1+f_m)^N / ((1+f_m)^N - 1); % monthly payment
    monthly_payment = g + h/12*p + z;
    str = ['Monthly payment when down payment is ' num2str(dp) ' is: ' num2str(monthly_payment)];
    disp(str);

    %% Compute mortgage ammortization schedule, equity, and loan balance

    loan_balance(1) = m;
    equity(1) = d;
    interest(1) = 0;

    for n = 1:N
        interest(n+1) = loan_balance(n)*(f_m);
        principle = g - interest(n+1);
        equity(n+1) = equity(n) + principle;
        loan_balance(n+1) = loan_balance(n) - principle;
    end

    if plot_mortgage && dp == 0.2
        figure;
        title('Owned value of home (ignoring appreciation) vs. Loan balance (20% down on 750k house)');
        hold on;
        plot((0:N)/12,equity,'-b');
        plot((0:N)/12,loan_balance,'-r');
        legend('Home equity', 'Loan balance');
        xlabel('Years');
        axis([0,30,0,1000000]);
        yticks([0,200000,400000,600000,800000,1000000]);
        yticklabels({'0','200k','400k','600k','800k','1m'});
        snapnow;
    end

    %% Compute value of ownership

    home_value(1) = p;
    owned_value(1) = d;
    fees(1) = 0;
    total_value(1) = d;
    for n = 1:N
        home_value(n+1) = home_value(n) * (1+(a-1)/12);
        owned_value(n+1) = (equity(n+1)/p) * home_value(n+1); % home value x amount owned
        fees(n+1) = fees(n) + z + interest(n+1);
        if mod(n,12) == 0
            fees(n+1) = fees(n+1) + h*home_value(n+1);
        end
        total_value(n+1) = owned_value(n+1) - fees(n+1);
    end

    if plot_ownership_value && dp == 0.2
        figure;
        title('Value analysis of owning a home (20% down on 750k house)');
        hold on;
        plot((0:N)/12,total_value,'-b');
        plot((0:N)/12,owned_value,'-g');
        plot((0:N)/12,fees,'-r');
        legend('Net value (Ownership)','Owned Home value', 'Taxes, fees and interest paid ');
        xlabel('Years');
        axis([0,30,0,3000000]);
        snapnow;
%         yticks([0,200000,400000,600000,800000,1000000]);
%         yticklabels({'0','200k','400k','600k','800k','1m'});
    end

    %% Compute value of renting / investing down payment

    trust_value(1) = d;
    rent(1) = r;
    rent_value(1) = d-r;

    for n = 1:N
        trust_value(n+1) = trust_value(n) * (1 + (b-1)/12) + (g+z+h/12*home_value(n+1)-r);
        rent(n+1) = rent(n) + r;
        rent_value(n+1) = trust_value(n+1)-rent(n+1);
    end

    if plot_renting_value && dp == 0.2
        figure;
        title('Value analysis of renting and keeping down payment in trust (20% down on 750k house)');
        hold on;
        plot((0:N)/12,rent_value,'-b');
        plot((0:N)/12,trust_value,'-g');
        plot((0:N)/12,rent,'-r');
        legend('Net value (Renting)','Trust value', 'Rent paid');
        xlabel('Years');
        axis([0,30,0,3000000]);
        snapnow;
%         yticks([0,200000,400000,600000,800000,1000000]);
%         yticklabels({'0','200k','400k','600k','800k','1m'});
    end
    %% Compare Owning vs. Renting/Investing
    value_diff{iter} = total_value-rent_value;
    iter = iter+1;
    
end

if plot_difference
    figure;
    title('Net Earnings When Owning Instead of Renting for different down payment percentages');
    xlabel('Years');
    hold on;
    for i = 1:numel(dps)
        values{i} = num2str(dps(i));
        plot((0:N)/12,value_diff{i});
    end
    values{numel(dps)+1} = 'Net equal';
    plot([0,30],[0,0],'--k');
    legend(values{:});
    snapnow;
end








