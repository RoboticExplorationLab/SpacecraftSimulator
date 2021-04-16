clear

run("parse_B.m")
run("parse_D.m")

% ecef_A = [-270073162	-429439276	385309511]/100;

lla = [37.3659105824291, -122.180314402498, 279.21]; % deg, deg, meters;

ecef_A = lla2ecef(lla);
%%
Rangeraw = RawRange;

idxs = find(~isnan(Rangeraw));

range = zeros(length(idxs),1);
dist = zeros(length(idxs),1);
time = zeros(length(idxs),1);
satsinfix = zeros(length(idxs),1);

for i = 1:length(idxs)
    
    % index with our range measurement 
    idx = idxs(i);
    
    % location of B
    ecef_B = B.og.ecef(:,idx-1)/100;
    satsinfix(i) = SVinFix1(idx-1);
    B.dist(i) = norm(ecef_ - ecef_B);
    
    % raw range measurement 
    B.range(i) = B.og.range(idx);
    
    % time 
    time(i) = TimeofWeek001s1(idx-1)/100;
end


figure
hold on 
yyaxis left
plot((time- time(1))/(60),dist)
ylabel('GPS Distance (meters)')

yyaxis right 
plot((time- time(1))/(60),satsinfix)
ylabel('Satellites in Fix')
xlabel('Time (minutes)')
hold off 
saveas(gcf,'satsinfix.png')


    
%% get stuff for D 
Rangeraw = RawRange1;
idxs = find(~isnan(Rangeraw));

range2 = zeros(length(idxs),1);
% dist = zeros(length(idxs),1);
time2 = zeros(length(idxs),1);

for i = 1:length(idxs)
    
    % index with our range measurement 
    idx = idxs(i);
    
    % raw range measurement 
    range2(i) = Rangeraw(idx);
    
    % time 
    time2(i) = TimeofWeek001s(idx-1)/100;
end


% time = (time - time(1));

gps_range = dist;
raw_range = range;
tow = time;
% save('max.mat','gps_range','raw_range','tow')

%%
% figure
% hold on
% plot(time,dist)
% xlabel('Time (s)')
% ylabel('GPS Distance (meters)')
% hold off 
% 
% figure
% hold on 
% plot(time,range)
% xlabel('Time (s)')
% ylabel('Raw Range Measurements')
% hold off 


%%
figure
hold on
yyaxis left
plot((time- time(1))/(60),dist)
ylabel('GPS Distance (meters)')
ylim([3860 3905])
hold off 
hold on 
yyaxis right
plot((time - time(1))/(60),range)
plot((time2- time(1))/(60),range2,'g')
ylabel('Raw Range Measurements')
% ylim([2.2e4 2.9e4] - 2e2)
ylim([1.085e5 1.096e5])
xlim([0 5.5])
xlabel('Time (minutes)')
hold off
saveas(gcf,'rangedata_secondsession_2.png')



%%
% yyaxis right 
% 
% xlabel('Time (s)')
% ylabel('GPS Distance (meters)')
% hold off 
% 
% figure
% hold on 
% plot(time,range)
% xlabel('Time (s)')
% ylabel('Raw Range Measurements')
% hold off 


% 
% 
% iter = 0
% dist = []
% rang
% for i = 1:length(Rangeraw)
%     r = Rangeraw(i);
%     if ~isnan(r)


%%

% lla = [37.402863, -122.165535, 136.39793] % deg, deg, meters;
figure
hold on 
plot(SVinFix)
plot(SVinFix1)
hold off 


%%
%% color plot for max 

% figure
% hold on 

% for i = 1:7
%     D{i} = dist;
%     SFidxs{i} = (satsinfix == i + 4);
%     Di = dist;
%     Di(~(satsinfix == i + 4)) = NaN;
%     D{i} = Di;
% %     D{i}
% end
% 
% 
% rgb1 = [29 38 113]/255;
%     rgb2 = [195 55 100]/255;
%     drgb = rgb2-rgb1;
%     
% figure
% hold on 
% for i = 1:7
%     plot((time- time(1))/(60),D{i},'Color',rgb1 + drgb*(i-1)/length(D))
% end
% hold off 

N = 1000 ;
x = (time- time(1))/(60);
y = dist;
z = satsinfix;
nx = linspace(x(1),x(end),N);
ny = spline(x,y,nx);
nz = spline(x,z,nx);

figure
hold on 
yyaxis right
plot((time - time(1))/(60),range)
plot((time2- time(1))/(60),range2)
ylim([1.085e5 1.096e5])
ylabel('Raw Range Data')


yyaxis left 
scatter(nx,ny,[],nz,'fill')
ylim([3860 3905])
% x = (time- time(1))/(60);
% y = dist;
% z = satsinfix;
% surf([x(:) x(:)], [y(:) y(:)], [z(:) z(:)], ...  % Reshape and replicate data
%      'FaceColor', 'none', ...    % Don't bother filling faces with color
%      'EdgeColor', 'interp', ...  % Use interpolated color for edges
%      'LineWidth', 2);            % Make a thicker line
c = colorbar;
c.Label.String = 'Satellites in Fix';
xlabel('Time (minutes)')
ylabel('GPS Distance (m)')


hold off 
% saveas(gcf,'colorplot.png')

%% least squares for affine fit 

t_B = (time2 - time(1))/(60);
t_D = (time - time(1))/(60);
lsx_B = range2 ;
lsx_D = range;
lsy = dist;

% figure
% hold on 
% plot(t_D,lsx_D)
% plot(t_B,lsx_B)
% hold off 


% trim after 2.5 
Dtrim = (t_D<2.5);
Btrim = (t_B<2.5);

t_D(Dtrim) = [];
t_B(Btrim) = [];
lsx_D(Dtrim) = [];
lsx_B(Btrim) = [];

% figure
% hold on 
% plot(t_D,lsx_D)
% plot(t_B,lsx_B)
% hold off 


% lsx_B = smoothdata(lsx_B);
% lsx_D = smoothdata(lsx_D);

A_B = [lsx_B ones(length(lsx_B),1)];
A_D = [lsx_D ones(length(lsx_D),1)];

x = (time- time(1))/(60);
y = dist;

y_B = spline(x,y,t_B);
y_D = spline(x,y,t_D);

xx = [A_B;A_D]\[y_B;y_D];
% xx = A_B\y_B;

figure
hold on 
plot(t_B,xx(1)*lsx_B + xx(2))
plot(t_D,xx(1)*lsx_D + xx(2))
plot(x,y)
hold off 

% %%
% figure
% hold on 
% plot(t_B,smoothdata(lsx_B))
% hold off 


