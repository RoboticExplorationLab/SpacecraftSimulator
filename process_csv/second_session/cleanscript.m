clear

run("parse_B.m")
run("parse_D.m")

lla = [37.3659105824291, -122.180314402498, 279.21]; % deg, deg, meters;

ecef_D = lla2ecef(lla)';

%% parse first data set 


idxs = find(~isnan(B.og.range));

B.range = zeros(length(idxs),1);
B.dist = zeros(length(idxs),1);
B.time = zeros(length(idxs),1);
B.SV = zeros(length(idxs),1);
B.FEE = zeros(length(idxs),1);

for i = 1:length(idxs)
    
    % index with our range measurement 
    idx = idxs(i);
    
    % location of B
    B.SV(i) = B.og.SV(idx-1);
    ecef_B = B.og.ecef(:,idx-1)/100;
    B.dist(i) = norm(ecef_D - ecef_B);
    
    % raw range measurement 
    B.range(i) = B.og.range(idx);
    
    % time 
    B.time(i) = B.og.time(idx-1)/100;
    
    % FEE
    B.FEE(i) = B.og.FEE(idx);
end
B.time = B.time - B.time(1);

%%
% figure
% hold on 
% plot(B.time,B.dist)
% hold off 
% 
% figure
% hold on 
% yyaxis left 
% plot(B.time,B.FEE)
% 
% yyaxis right 
% plot(B.time,B.range)
% hold off 

figure
hold on 
subplot(2,1,1)
plot(B.time,B.dist)
ylabel('GPS Distance (m)')
xlabel('Time (min)')

subplot(2,1,2)
plot(B.time,B.FEE)
ylabel('Frequency Offset')
xlabel('Time (min)')
saveas(gcf,'dist_fee.png')
hold off

%% least squares time
% trim after 2.5 
Btrim = (B.time<(2.5*60));

B.trim.time = B.time;
B.trim.time(Btrim) = [];
B.trim.range = B.range;
B.trim.range(Btrim) = [];
B.trim.dist = B.dist;
B.trim.dist(Btrim) = [];
B.trim.FEE = B.FEE;
B.trim.FEE(Btrim) = [];

ls_sol = [B.trim.range ones(length(B.trim.range),1)]\B.trim.dist;
ls_sol_fee = [B.trim.range ones(length(B.trim.range),1) B.trim.FEE]\B.trim.dist;

figure
hold on 
title('Linear Regression (No Pre-Smoothing)')
plot(B.trim.time/60,B.trim.dist)
plot(B.trim.time/60,ls_sol(1)*B.trim.range + ls_sol(2))
plot(B.trim.time/60,ls_sol_fee(1)*B.trim.range + ls_sol_fee(2) + ls_sol_fee(3)*B.trim.FEE)
ylim([3860 3905])
ylabel('Distance (m)')
xlabel('Time (min)')
legend('GPS Distance (m)','Transformed Range (No FEE)','Transformed Range (w/ FEE)')
saveas(gcf,'LS_nosmooth.png')
hold off 

% %% Huber version 
% Btrim = (B.time<(2.5*60));
% 
% B.trim.time = B.time;
% B.trim.time(Btrim) = [];
% B.trim.range = B.range;
% B.trim.range(Btrim) = [];
% B.trim.dist = B.dist;
% B.trim.dist(Btrim) = [];
% B.trim.FEE = B.FEE;
% B.trim.FEE(Btrim) = [];
% 
% ls_sol = robustfit(B.trim.range,B.trim.dist);
% % ls_sol_fee = robustfit(B.trim.range ones(length(B.trim.range),1) B.trim.FEE]\B.trim.dist;
% 
% figure
% hold on 
% title('Linear Regression (No Pre-Smoothing)')
% plot(B.trim.time/60,B.trim.dist)
% plot(B.trim.time/60,ls_sol(2)*B.trim.range + ls_sol(1))
% % plot(B.trim.time/60,ls_sol_fee(1)*B.trim.range + ls_sol_fee(2) + ls_sol_fee(3)*B.trim.FEE)
% % ylim([3860 3905])
% ylabel('Distance (m)')
% xlabel('Time (min)')
% legend('GPS Distance (m)','Transformed Range (No FEE)','Transformed Range (w/ FEE)')
% % saveas(gcf,'LS_nosmooth.png')
% hold off 

%% smoothed version 
% trim after 2.5 
Btrim = (B.time<(2.5*60));

B.trim.time = B.time;
B.trim.time(Btrim) = [];
B.trim.range = B.range;
B.trim.range(Btrim) = [];
B.trim.dist = B.dist;
B.trim.dist(Btrim) = [];
B.trim.FEE = B.FEE;
B.trim.FEE(Btrim) = [];

B.trim.range = smoothdata(B.trim.range,'movmean',10);
ls_sol = [(B.trim.range) ones(length(B.trim.range),1)]\B.trim.dist;
% ls_sol = [smoothdata(B.trim.range) ones(length(B.trim.range),1)]\B.trim.dist;
% ls_sol_fee = [smoothdata(B.trim.range) ones(length(B.trim.range),1) smoothdata(B.trim.FEE)]\B.trim.dist;
e = [(B.trim.range) ones(length(B.trim.range),1)]*ls_sol - B.trim.dist;
RMS = sqrt(mean( e .* e ))
figure
hold on 
title('Linear Regression (With Pre-Smoothing)')
plot(B.trim.time/60,B.trim.dist)
plot(B.trim.time/60,ls_sol(1)*B.trim.range + ls_sol(2))
% plot(B.trim.time/60,ls_sol_fee(1)*B.trim.range + ls_sol_fee(2) + ls_sol_fee(3)*B.trim.FEE)
ylim([3860 3905])
ylabel('Distance (m)')
xlabel('Time (min)')
legend('GPS Distance (m)','Transformed Range (No FEE)','Transformed Range (w/ FEE)')
saveas(gcf,'LS_smooth.png')
hold off 

%% new nice plot 
close all
figure
hold on 
% pltsigma = sqrt(2)*2.5*2.0789; % meters
pltsigma = sqrt(2)*2.5*2.89; % meters
% title('Ranging Calibration (B + E)')
title(['Ranging Calibration (B to D) ' '   (RMS = ' num2str(RMS) 'm)'],'FontSize',24)
% newt = (E.interp.time - E.interp.time(1))/60;
%     plot(newt,Erange/1000,'linewidth',2)
    [hl, ~]  = boundedline(B.trim.time/60,B.trim.dist/1000,pltsigma/1000);
    hl.LineWidth = 2;
    plot(B.trim.time/60,(ls_sol(1)*B.trim.range + ls_sol(2))/1000,'linewidth',2,'Color',[0.8500 0.3250 0.0980])
%     plot(newt,smoothdata(E_coeffs(1)*E.interp.range + E_coeffs(2),'movmedian',smoothwindow)/1000,'linewidth',2,'Color',[0.8500 0.3250 0.0980])
    ylabel('Range (km)','FontSize',14)
xlabel('Time (minutes)','FontSize',14)
% xlim([newt(1),newt(end)])
grid on 
legend('GPS Accuracy 3\sigma','GPS Distance','Calibrated Range','FontSize',18,'Location','best')
ax = gca(); 
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
saveas(gcf,'first_test.svg')

hold off 
% close all

%%

% ee = [smoothdata(B.trim.range) ones(length(B.trim.range),1)]*ls_sol - B.trim.dist;
ee = [(B.trim.range) ones(length(B.trim.range),1)]*ls_sol - B.trim.dist;
sqrt(mean( ee .* ee ))

figure
hold on 
title('position errors after affine fit')
plot(B.trim.time/60,ee)
xlabel('Time (min)')
ylabel('Error (m)')
saveas(gcf,'error.png')

%% cross correlation 

r = xcorr(B.trim.FEE,B.trim.dist);
r = xcorr(B.trim.FEE,B.trim.range);
hold off 

