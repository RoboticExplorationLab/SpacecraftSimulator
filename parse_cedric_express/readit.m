%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/kevintracy/Downloads/vr3x-balloon-tracker-data.xlsx
%    Worksheet: vr3x-balloon-tracker-data
%
% Auto-generated by MATLAB on 15-Mar-2021 14:13:40

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "vr3x-balloon-tracker-data";
opts.DataRange = "C2:G18843";

% Specify column names and types
opts.VariableNames = ["DateTimeCST", "User", "Latitude", "Longitude", "Altitudem"];
opts.SelectedVariableNames = ["DateTimeCST", "User", "Latitude", "Longitude", "Altitudem"];
opts.VariableTypes = ["datetime", "string", "double", "double", "double"];
opts = setvaropts(opts, 1, "InputFormat", "");
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");

% Import the data
vr3xballoontrackerdata = readtable("/Users/kevintracy/Downloads/vr3x-balloon-tracker-data.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts

DATA = vr3xballoontrackerdata;

%% separate 

adam.data = DATA(1:2821,:);

alex.data = DATA(2822:4937,:);

anh.data = DATA(4938:6273,:);

baloon.data = DATA(6274:16112,:);

cedric.data = DATA(16113:18288,:);

kevin.data = DATA(18289:end,:);

C = {adam,alex,anh,baloon,cedric,kevin};

