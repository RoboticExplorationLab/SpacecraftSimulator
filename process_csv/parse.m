clear

run("parseB.m")


% ecef_A = [-270073162	-429439276	385309511]/100;

lla = [37.402863, -122.165535, 136.39793]; % deg, deg, meters;

ecef_A = lla2ecef(lla);
%%


idxs = find(~isnan(Rangeraw));

range = zeros(length(idxs),1);
dist = zeros(length(idxs),1);
time = zeros(length(idxs),1);

for i = 1:length(idxs)
    
    % index with our range measurement 
    idx = idxs(i);
    
    % location of B
    ecef_B = [ECEFX001m(idx-1),ECEFY001m(idx-1),ECEFZ001m(idx-1)]/100;
    dist(i) = norm(ecef_A - ecef_B);
    
    % raw range measurement 
    range(i) = Rangeraw(idx);
    
    % time 
    time(i) = TimeofWeek001s(idx-1)/100;
end


time = (time - time(1));

gps_range = dist
raw_range = range
tow = time
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
plot(time,dist)
ylabel('GPS Distance (meters)')
hold off 
hold on 
yyaxis right
plot(time,range)
ylabel('Raw Range Measurements')
ylim([2.2e4 2.9e4] - 2e2)
xlabel('Time (minutes)')
hold off
saveas(gcf,'rangedata.png')



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

