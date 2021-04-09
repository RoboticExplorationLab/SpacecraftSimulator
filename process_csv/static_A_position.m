clear



run("parseA.m")

%%


figure
hold on 
% yyaxis left
plot(trimnan(ECEFX001m1)/100)
ylabel('ECEF X (m)')

% yyaxis right 
% plot(trimnan(SVinFix1))
% ylabel('Sats in Fix')
% plot(trimnan(ECEFY001m1))
% plot(trimnan(ECEFZ001m1))
hold off 


min(trimnan(ECEFX001m1)/100)
max(trimnan(ECEFX001m1)/100)




r1 = [ECEFX001m1(400), ECEFY001m1(400), ECEFZ001m1(400)]/100;


function Aout =  trimnan(A)
Aout = A ;
Aout(A == 0) = NaN;
end