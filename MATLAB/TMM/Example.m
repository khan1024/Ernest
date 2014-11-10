%% Initialization

%House keeping
clc;
clear all;

%Sets beginning hours and minutes
Time.Hours =9;
Time.Minutes =28;

%Sets the day of year (1 to 365)
Time.DayOfYear = 82;

%Time zone (see "help SolarRadiation")
Time.StandardMeridian = 8;

%Number of days to run for
Time.RunLength = 90;

%Latitude and longitude
Location.Latitude = 52;
Location.Longitude = -115;

%Tilt angle of the array (use 0 degrees for default)
Location.TiltAngle = 0;

%Ground reflectance, for diffuse light calculations (default is 0.2)
Location.GroundReflectance = 0.2;

%Put "1" for no clouds (AM1.5), "0" for full obscuration
Location.ClearnessFactor = 1;


%% Runs the script and plots results

%Performs the solar power and angle calculation
Results = SolarRadiation(Time, Location);

AvgW = (sum(Results.DirectRadiation)+sum(Results.DiffuseRadation))/(1440*Time.RunLength)
AvgkWhrPerDay = AvgW * 24 / 1000


%Plots solar radation
subplot(3,1,1);
plot(Results.Time/60, Results.DirectRadiation);
xlabel('Time (hrs)');
ylabel('Direct radiation (W/m^2)');


%Plots the solar azimuth angle
subplot(3,1,2);
plot(Results.Time/60, Results.SunAzimuthAngle);
xlabel('Time (hrs)');
ylabel('Azimuth angle (degrees)');


%Plots the solar elevation angle
subplot(3,1,3);
plot(Results.Time/60, Results.SunElevationAngle);
xlabel('Time (hrs)');
ylabel('elevation angle (degrees)');
