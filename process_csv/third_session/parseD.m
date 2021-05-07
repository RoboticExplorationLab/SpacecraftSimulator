function [D] = parseD(t_vec)
opts = delimitedTextImportOptions("NumVariables", 21);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["FixMode", "SVinFix", "WeekNum", "TimeofWeek001s", "Latitude1E7deg", "Longitude1E7deg", "HeightAboveEllipsoid001m", "HeightAboveSeaLevel001m", "GeometricDilutionofPrecision001", "PositionDilutionofPrecision001", "HorizontalDilutionofPrecision001", "VerticalDilutionofPrecision001", "TimeDilutionofPrecision001", "ECEFX001m", "ECEFY001m", "ECEFZ001m", "ECEFVX001ms", "ECEFVY001ms", "ECEFVZ001ms", "RawRange", "FrequencyOffset"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("/Users/kevintracy/devel/SpacecraftSimulator/process_csv/third_session/4D_range_w_gps.csv", opts);

%% Convert to output type
FixMode = tbl.FixMode;
SVinFix = tbl.SVinFix;
WeekNum = tbl.WeekNum;
TimeofWeek001s = tbl.TimeofWeek001s;
Latitude1E7deg = tbl.Latitude1E7deg;
Longitude1E7deg = tbl.Longitude1E7deg;
HeightAboveEllipsoid001m = tbl.HeightAboveEllipsoid001m;
HeightAboveSeaLevel001m = tbl.HeightAboveSeaLevel001m;
GeometricDilutionofPrecision001 = tbl.GeometricDilutionofPrecision001;
PositionDilutionofPrecision001 = tbl.PositionDilutionofPrecision001;
HorizontalDilutionofPrecision001 = tbl.HorizontalDilutionofPrecision001;
VerticalDilutionofPrecision001 = tbl.VerticalDilutionofPrecision001;
TimeDilutionofPrecision001 = tbl.TimeDilutionofPrecision001;
ECEFX001m = tbl.ECEFX001m;
ECEFY001m = tbl.ECEFY001m;
ECEFZ001m = tbl.ECEFZ001m;
ECEFVX001ms = tbl.ECEFVX001ms;
ECEFVY001ms = tbl.ECEFVY001ms;
ECEFVZ001ms = tbl.ECEFVZ001ms;
RawRange = tbl.RawRange;
FrequencyOffset = tbl.FrequencyOffset;


for i = 1:length(ECEFX001m)
    B.og.ecef(:,i) = [ECEFX001m(i);ECEFY001m(i);ECEFZ001m(i)];
end
B.og.range= RawRange;
B.og.FEE = FrequencyOffset;
B.og.SV = SVinFix;
B.og.time = TimeofWeek001s;


% B.nu.ecef = B.og.ecef;
% B.nu.ecef(:,isnan(B.og.ecef(1,:))) = [];
% B.nu.time = B.og.time;
% B.nu.time(isnan(B.og.ecef(1,:))) = [];
% B.nu.SVinFix = SVinFix;
% B.nu.SVinFix(isnan(B.og.ecef(1,:))) = [];
% r0 = B.nu.ecef(:,1);
% norm(r0)
% for i = 1:length(B.nu.time)
%     err(i) = norm(r0 - B.nu.ecef(:,i))/100;
% end

% figure
% hold on 
% title('Static GPS Deviations','FontSize',24)
% plot((B.nu.time - B.nu.time(1))/(60*100),err,'linewidth',2)
% grid on
% xlabel('Time (minutes)')
% ylabel('Distance From First Measurement (m)')
% ax = gca(); 
% ax.XAxis.FontSize = 16;
% ax.YAxis.FontSize = 16;
% % saveas(gcf,'gps_errors.svg')
% hold off 
% 
% figure
% hold on 
% title('Satellites in Fix for Static Deviations','FontSize',24)
% plot((B.nu.time - B.nu.time(1))/(60*100),B.nu.SVinFix,'linewidth',2)
% grid on
% xlabel('Time (minutes)')
% ylabel('Distance From First Measurement (m)')
% ax = gca(); 
% ax.XAxis.FontSize = 16;
% ax.YAxis.FontSize = 16;
% % saveas(gcf,'gps_errors.svg')
% hold off 


% 
% xx = (B.nu.time - B.nu.time(1))/(60*100);
% figure
% set(gca,'FontSize',16)
% hold on 
% ax = gca(); 
% ax.XAxis.FontSize = 16;
% ax.YAxis.FontSize = 16;
% sg = sgtitle('Static GPS Performance');
% sg.FontSize = 24
% title('test')
% subplot(2,1,1)
% % set(gca,'FontSize',16)
% plot((B.nu.time - B.nu.time(1))/(60*100),err,'linewidth',2)
% grid on
% xlabel('Time (minutes)')
% ylabel('\Delta Position (m)')
% xlim([xx(1),xx(end)])
% % ax = gca(); 
% % ax.XAxis.FontSize = 16;
% % ax.YAxis.FontSize = 16;
% ax = gca(); 
% ax.XAxis.FontSize = 16;
% ax.YAxis.FontSize = 16;
% 
% subplot(2,1,2)
% hold on 
% plot((B.nu.time - B.nu.time(1))/(60*100),B.nu.SVinFix,'linewidth',2)
% grid on
% xlabel('Time (minutes)')
% ylabel('Satellites in Fix')
% ax = gca(); 
% ax.XAxis.FontSize = 16;
% ax.YAxis.FontSize = 16;
% xlim([xx(1),xx(end)])
% yticks([8 9 10 11 12])
% % set(findall(gcf,'-property','FontSize'),'FontSize',16)
% % fh = findall(0,'Type','Figure');
% % set( findall(fh, '-property', 'fontsize'), 'fontsize', 14)
% hold off
% saveas(gcf,'gps_errors2.svg')
% % ax = gca(); 
% % ax.XAxis.FontSize = 16;
% % ax.YAxis.FontSize = 16;












mat = [B.og.time B.og.range];
% remove NaN's
mat(isnan(mat(:,2)),:) = [];
% remove duplicate values 
mat((diff(mat(:,1)) == 0 ),:) = [];

% t_vec = linspace(5.04e5,5.05e5,1000);

% convert time to seconds 
mat(:,1) = mat(:,1)/100;

% create interp objects
B.interp.time = t_vec;
% B.interp.range = vec(spline(mat(:,1),mat(:,2),t_vec));
B.interp.range = vec(interp1(mat(:,1),mat(:,2),t_vec));
% get name right
D = B;
end