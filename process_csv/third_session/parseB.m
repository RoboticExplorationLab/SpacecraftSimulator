function [B] = parseB()
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
tbl = readtable("/Users/kevintracy/devel/SpacecraftSimulator/process_csv/third_session/4B_range_w_gps.csv", opts);
tbl = tbl(13:end,:); % trim the first values
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

% get name right
% E = B;
end