clear

B = parseB();

t_vec = linspace(5.04e5,5.05e5,1000);
% tf
%% create interp objects
ecef_mat = [B.og.time B.og.ecef'];

% remove NaN's
idxs = find(isnan(B.og.ecef(1,:)));
ecef_mat(idxs,:) = [];

% remove duplicate values 
ecef_mat((diff(ecef_mat(:,1)) == 0 ),:) = [];

% put in propper structs
B.ecef.time = ecef_mat(:,1)/100;
B.ecef.data = ecef_mat(:,2:4)'/100;

% figure
% hold on 
% plot(B.ecef.time,B.ecef.data(1,:))
% hold off 

B.interp.time = t_vec;
B.interp.ecef(1,:) = spline(B.ecef.time,B.ecef.data(1,:),t_vec);
B.interp.ecef(2,:) = spline(B.ecef.time,B.ecef.data(2,:),t_vec);
B.interp.ecef(3,:) = spline(B.ecef.time,B.ecef.data(3,:),t_vec);

ecef_A = [-270137265,	-429420935,	385286916]'/100;
ecef_C = [-270668368,	-429743088,	384652272]'/100;
ecef_D = [-270831340,	-429305807,	385001248]'/100;
ecef_E = [-270322037,	-429590411,	384989518]'/100;

Erange = zeros(size(B.interp.ecef,2),1);
Crange = zeros(size(B.interp.ecef,2),1);
Drange = zeros(size(B.interp.ecef,2),1);
Arange = zeros(size(B.interp.ecef,2),1);
for i = 1:size(B.interp.ecef,2)
    r = B.interp.ecef(:,i);
    Erange(i) = norm(vec(ecef_E) - vec(r));
    Crange(i) = norm(ecef_C - r);
    Drange(i) = norm(ecef_D - r);
    Arange(i) = norm(ecef_A - r);
end

% figure
% hold on 
% plot(t_vec,Erange)
% hold off 

% E = parseE();

% figure
% hold on 
% plot(E.interp.time,E.interp.range)
% hold off 

%% least squares 

E = parseE();
performLS(Erange,E,'Ranging Calibration (B + E)','BE.png')

C = parseC();
performLS(Crange,C,'Ranging Calibration (B + C)','BC.png')

D = parseD();
performLS(Drange,D,'Ranging Calibration (B + D)','BD.png')

A = parseA();
performLS(Arange,A,'Ranging Calibration (B + A)','BA.png')





%% 
function [] = performLS(Erange,E,titlename,filename)
E_coeffs = [E.interp.range ones(1000,1) ]\Erange;
e = [E.interp.range ones(1000,1) ]*E_coeffs - Erange;
display(filename)
RMS  = sqrt( mean ( e .* e ) )
close all
figure
hold on 
% title('Ranging Calibration (B + E)')
title(titlename)
plot(E.interp.time,Erange)
plot(E.interp.time,E_coeffs(1)*E.interp.range + E_coeffs(2))
xlabel('GPS TOW (s)')
ylabel('Range (m)')
legend('GPS Distance','Calibrated Range Measurements')
saveas(gcf,filename)
% saveas(gcf,'test.png','epsc')
hold off 
close all
end

%% 
% %%
% figure
% hold on 
% plot(B.ecef.time)
% hold off 
% 
% B.ecef.time((diff(B.ecef.time) == 0 )) = [];
% figure
% hold on 
% plot(B.ecef.time)
% hold off 
% 
% 
% %%
% % idxs = find(~isnan(B.og.range));
% 
% % B.range = zeros(length(idxs),1);
% % B.dist = zeros(length(idxs),1);
% % B.time = zeros(length(idxs),1);
% % B.SV = zeros(length(idxs),1);
% % B.FEE = zeros(length(idxs),1);
% % B.ecef = zeros(3,length(idxs));
% % 
% % for i = 1:length(idxs)
% %     
% %     % index with our range measurement 
% %     idx = idxs(i);
% %     
% %     % location of B
% %     B.SV(i) = B.og.SV(idx-1);
% %     ecef_B = B.og.ecef(:,idx-1)/100;
% %     B.ecef(:,i) = ecef_B;
% % %     B.dist(i) = norm(ecef_D - ecef_B);
% %     
% %     % raw range measurement 
% %     B.range(i) = B.og.range(idx);
% %     
% %     % time 
% %     B.time(i) = B.og.time(idx-1)/100;
% %     
% %     % FEE
% %     B.FEE(i) = B.og.FEE(idx);
% % end
% % B.time = B.time - B.time(1);
% 
% %%
% % figure
% % hold on 
% % plot(B.range)
% % plot(smoothdata(B.range,'movmedian',7))
% % hold off 
% 
% % figure
% % hold on 
% % plot(B.ecef(1,:))
% % hold off 
% 
% %%
% % figure
% % hold on 
% % plot(B.time,B.dist)
% % hold off 
% % 
% % figure
% % hold on 
% % yyaxis left 
% % plot(B.time,B.FEE)
% % 
% % yyaxis right 
% % plot(B.time,B.range)
% % hold off 
% 
% figure
% hold on 
% subplot(2,1,1)
% plot(B.time,B.dist)
% ylabel('GPS Distance (m)')
% % xlabel('Time (min)')
% 
% subplot(2,1,2)
% plot(B.time,B.FEE)
% ylabel('Frequency Offset')
% xlabel('Time (min)')
% % saveas(gcf,'dist_fee.png')
% hold off
% 
% %% least squares time
% % trim after 2.5 
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
% ls_sol = [B.trim.range ones(length(B.trim.range),1)]\B.trim.dist;
% ls_sol_fee = [B.trim.range ones(length(B.trim.range),1) B.trim.FEE]\B.trim.dist;
% 
% figure
% hold on 
% title('Linear Regression (No Pre-Smoothing)')
% plot(B.trim.time/60,B.trim.dist)
% plot(B.trim.time/60,ls_sol(1)*B.trim.range + ls_sol(2))
% plot(B.trim.time/60,ls_sol_fee(1)*B.trim.range + ls_sol_fee(2) + ls_sol_fee(3)*B.trim.FEE)
% ylim([3860 3905])
% ylabel('Distance (m)')
% xlabel('Time (min)')
% legend('GPS Distance (m)','Transformed Range (No FEE)','Transformed Range (w/ FEE)')
% saveas(gcf,'LS_nosmooth.png')
% hold off 
% 
% % %% Huber version 
% % Btrim = (B.time<(2.5*60));
% % 
% % B.trim.time = B.time;
% % B.trim.time(Btrim) = [];
% % B.trim.range = B.range;
% % B.trim.range(Btrim) = [];
% % B.trim.dist = B.dist;
% % B.trim.dist(Btrim) = [];
% % B.trim.FEE = B.FEE;
% % B.trim.FEE(Btrim) = [];
% % 
% % ls_sol = robustfit(B.trim.range,B.trim.dist);
% % % ls_sol_fee = robustfit(B.trim.range ones(length(B.trim.range),1) B.trim.FEE]\B.trim.dist;
% % 
% % figure
% % hold on 
% % title('Linear Regression (No Pre-Smoothing)')
% % plot(B.trim.time/60,B.trim.dist)
% % plot(B.trim.time/60,ls_sol(2)*B.trim.range + ls_sol(1))
% % % plot(B.trim.time/60,ls_sol_fee(1)*B.trim.range + ls_sol_fee(2) + ls_sol_fee(3)*B.trim.FEE)
% % % ylim([3860 3905])
% % ylabel('Distance (m)')
% % xlabel('Time (min)')
% % legend('GPS Distance (m)','Transformed Range (No FEE)','Transformed Range (w/ FEE)')
% % % saveas(gcf,'LS_nosmooth.png')
% % hold off 
% 
% %% smoothed version 
% % trim after 2.5 
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
% ls_sol = [smoothdata(B.trim.range) ones(length(B.trim.range),1)]\B.trim.dist;
% ls_sol_fee = [smoothdata(B.trim.range) ones(length(B.trim.range),1) smoothdata(B.trim.FEE)]\B.trim.dist;
% 
% figure
% hold on 
% title('Linear Regression (With Pre-Smoothing)')
% plot(B.trim.time/60,B.trim.dist)
% plot(B.trim.time/60,ls_sol(1)*B.trim.range + ls_sol(2))
% plot(B.trim.time/60,ls_sol_fee(1)*B.trim.range + ls_sol_fee(2) + ls_sol_fee(3)*B.trim.FEE)
% ylim([3860 3905])
% ylabel('Distance (m)')
% xlabel('Time (min)')
% legend('GPS Distance (m)','Transformed Range (No FEE)','Transformed Range (w/ FEE)')
% saveas(gcf,'LS_smooth.png')
% hold off 
% 
% %%
% 
% % ee = [smoothdata(B.trim.range) ones(length(B.trim.range),1)]*ls_sol - B.trim.dist;
% ee = [(B.trim.range) ones(length(B.trim.range),1)]*ls_sol - B.trim.dist;
% sqrt(mean( ee .* ee ))
% 
% figure
% hold on 
% title('position errors after affine fit')
% plot(B.trim.time/60,ee)
% xlabel('Time (min)')
% ylabel('Error (m)')
% saveas(gcf,'error.png')
% 
% %% cross correlation 
% 
% r = xcorr(B.trim.FEE,B.trim.dist);
% r = xcorr(B.trim.FEE,B.trim.range);
% hold off 
% 
% 
% function time =  calibrate_time(time,t0)
% 
%     time = time - t0 ;
%     time = time/100;
% %     time = time/60 ;
% end