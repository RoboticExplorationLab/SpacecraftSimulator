clear

run("import_full.m")

%%
adam.data = DATA(1:2821,:);

alex.data = DATA(2822:4937,:);

anh.data = DATA(4938:6273,:);

balloon.data = DATA(6275:16112,:);

cedric.data = DATA(16113:18288,:);

kevin.data = DATA(18289:end,:);

C = {adam,alex,anh,balloon,cedric,kevin};

names = ["adam","alex","anh","baloon","cedric","kevin"];
t_balloon0 = balloon.data(1,1);

for i = 1:length(C)
    C{i}.t = ((C{i}.data(:,1) - t_balloon0)/1000)/60;
end



%% trim values after landing 
t_land = ((1615571828290 - t_balloon0)/1000)/60;

for i = 1:length(C)
    for j = size(C{i}.data,1):(-1):1
        if C{i}.t(j)>t_land
            C{i}.t(j) = [];
            C{i}.data(j,:) = [];
%             C{i}.data(j,5) = [];
%             C{i}.data(j,5) = [];
        end
    end
end
            
%%
figure
hold on 
for i = 1:length(C)
%     C{i}.t = C{i}.data(:,1);
%     plot(C{i}.t,C{i}.data(:,5)
long = C{i}.data(:,6);
lat = C{i}.data(:,5);

size(long)
size(lat)
plot(long,lat)
% plot(1:20,randn(20,1))
% plot(C{i}.data(:,6),C{i}.data(:,5))
end
legend(names)
hold off 


for i = 1:length(C)
    D = C{i};
    t = D.t;
    long = C{i}.data(:,6);
    lat = C{i}.data(:,5);
    
figure
hold on 
plot(t,lat,'o')
% plot(t,long)
title(names(i))
hold off 
end


%%

figure
hold on 
for i = 1:length(C)
%     C{i}.t = C{i}.data(:,1);
%     plot(C{i}.t,C{i}.data(:,5)
t = C{i}.t;
long = C{i}.data(:,6);
lat = C{i}.data(:,5);

size(long)
size(lat)
plot(t,lat,'o')
% plot(1:20,randn(20,1))
% plot(C{i}.data(:,6),C{i}.data(:,5))
end
legend(names)
hold off 



%% now the annoying work 

% adam 
% 0 - 98.2 (43.8424,-95.1555)
% 192 - 288 (43.6517,-93.7230)


% alex 
% 0 - 71 (44.0076799,-95.8058818)
% 151 - 288  (43.9640347,-94.4233140)


% anh 

% balloon

% cedric 

% kevin 
% 0 - 62 (43.646739, -95.993062)
% 159 - 288 (43.672691, -94.207618)


%% create new time vectors for ground stations 

t = C{4}.t;

alex_vec = NaN*zeros(3,length(t),1);
kevin_vec = NaN*zeros(3,length(t),1);
adam_vec = NaN*zeros(3,length(t),1);

for i = 1:length(t)
    
    tc = t(i);
    
    if tc <62
        kevin_vec(:,i) = [43.646739, -95.993062, 450];
    end
    if tc > 159
        kevin_vec(:,i) = [43.672691, -94.207618, 450];
    end
    
    if tc <71
        alex_vec(:,i) = [44.0076799,-95.8058818, 450];
    end
    if tc > 151
        alex_vec(:,i) = [43.9640347,-94.4233140, 450];
    end
    
    if tc <98.2
        adam_vec(:,i) = [43.8424,-95.1555, 450];
    end
    if tc > 192
        adam_vec(:,i) = [43.6517,-93.7230, 450];
    end
    
end


%% now get the distances

kevin_range = zeros(length(t),1);
adam_range = zeros(length(t),1);
alex_range = zeros(length(t),1);
cedric_range = zeros(length(t),1);

for i = 1:length(t)
    r_balloon = lla2ecef(C{4}.data(i,5:7));
    
    r_kevin = lla2ecef(kevin_vec(:,i)');
    r_alex = lla2ecef(alex_vec(:,i)');
    r_adam = lla2ecef(adam_vec(:,i)');
    
    
    r_cedric = lla2ecef(C{5}.data(1,5:7));
    
    kevin_range(i) = norm(r_balloon - r_kevin)/1000;
    alex_range(i) = norm(r_balloon - r_alex)/1000;
    adam_range(i) = norm(r_balloon - r_adam)/1000;
    cedric_range(i) = norm(r_balloon - r_cedric)/1000;
    
end

figure
hold on 
plot(t,kevin_range,'linewidth',2)
plot(t,adam_range,'linewidth',2)
plot(t,alex_range,'linewidth',2)
plot(t,cedric_range,'linewidth',2)
ylabel('Range (km)')
xlabel('Flight Time (min)')
grid on 
legend('Unit A','Unit C','Unit D','Unit E','location','northwest')
hold off 
% saveas(gcf,'rangeplot.png')
saveas(gcf,'rangeplot.eps','epsc')

%%
OUTDATA = [t kevin_range alex_range adam_range cedric_range];


%%
balloon_alt = C{4}.data(:,7);
figure
hold on 
plot([0;t + 8],[431.9;balloon_alt]/1000,'linewidth',2)
grid on 
ylabel('Altitude (km)')
xlabel('Flight Time (min)')
hold off
saveas(gcf,'altplot.png')
% saveas(gcf,'altplot.eps','epsc')