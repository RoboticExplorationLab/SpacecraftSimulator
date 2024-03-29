%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/kevintracy/devel/SpacecraftSimulator/process_csv/4A_range_w_gps.txt
%
% Auto-generated by MATLAB on 08-Apr-2021 12:53:14

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 20);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["FixMode1", "SVinFix1", "WeekNum1", "TimeofWeek001s1", "Latitude1E7deg1", "Longitude1E7deg1", "HeightAboveEllipsoid001m1", "HeightAboveSeaLevel001m1", "GeometricDilutionofPrecision1", "PositionDilutionofPrecision1", "HorizontalDilutionofPrecision1", "VerticalDilutionofPrecision1", "TimeDilutionofPrecision1", "ECEFX001m1", "ECEFY001m1", "ECEFZ001m1", "ECEFVX001ms1", "ECEFVY001ms1", "ECEFVZ001ms1", "Rangeraw1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
opts = setvaropts(opts, 20, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 20, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("/Users/kevintracy/devel/SpacecraftSimulator/process_csv/4A_range_w_gps.txt", opts);

%% Convert to output type
FixMode1 = tbl.FixMode1;
SVinFix1 = tbl.SVinFix1;
WeekNum1 = tbl.WeekNum1;
TimeofWeek001s1 = tbl.TimeofWeek001s1;
Latitude1E7deg1 = tbl.Latitude1E7deg1;
Longitude1E7deg1 = tbl.Longitude1E7deg1;
HeightAboveEllipsoid001m1 = tbl.HeightAboveEllipsoid001m1;
HeightAboveSeaLevel001m1 = tbl.HeightAboveSeaLevel001m1;
GeometricDilutionofPrecision1 = tbl.GeometricDilutionofPrecision1;
PositionDilutionofPrecision1 = tbl.PositionDilutionofPrecision1;
HorizontalDilutionofPrecision1 = tbl.HorizontalDilutionofPrecision1;
VerticalDilutionofPrecision1 = tbl.VerticalDilutionofPrecision1;
TimeDilutionofPrecision1 = tbl.TimeDilutionofPrecision1;
ECEFX001m1 = tbl.ECEFX001m1;
ECEFY001m1 = tbl.ECEFY001m1;
ECEFZ001m1 = tbl.ECEFZ001m1;
ECEFVX001ms1 = tbl.ECEFVX001ms1;
ECEFVY001ms1 = tbl.ECEFVY001ms1;
ECEFVZ001ms1 = tbl.ECEFVZ001ms1;
Rangeraw1 = tbl.Rangeraw1;

%% Clear temporary variables
clear opts tbl