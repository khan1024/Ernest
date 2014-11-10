function Results = SolarRadiation(Time, Location)
% function Results = SolarRadiation(Time, Location)
% 
% Usage: Projects solar radiation and angles given a location and time
% 
% Inputs:
% 
%     Time.Hours -> starting hour     
%     Time.Minutes -> starting minute
%     Time.DayOfYear -> starting day of year
%     Time.StandardMeridian -> time zone (see lookup table below)
%     Time.RunLength -> run length (days)
%     
%     Location.Latitude -> latitude (decimal degrees)
%     Location.Longitude -> longitude (decimal degrees)
%     Location.TiltAngle -> Earth's tilt angle (degrees)
%     Location.GroundReflectance -> ground reflectance (typically 0.2) 
%     Location.ClearnessFactor -> sky clearness factor (typically 1)
%
%     Standard Meridian Lookup Table:
% 
%     0 = 0 Greenwich Mean T. (0 deg); London
%     1 = -1 West African T. (15 deg); Azoren
%     2 = -2 Azores T. (30 deg); Fernando de Noronha
%     3 = -3 Eastern Brazilian T. (45 deg); Rio de Janeiro
%     4 = -4 Atlantic S.T. (60 deg); La Paz
%     5 = -5 Eastern S.T. (75 deg); New York
%     6 = -6 Central (90 deg); Chicago
%     7 = -7 Mountain (105 deg); Denver
%     8 = -8 Pacific (120 deg); Los Angeles
%     9 = -9  Alaska S.T. (135 deg); Anchorage
%     10 = -10  Alaska-Hawaii S.T. (150 deg); Honolulu
%     11 = -11 Samoa S.T. (165 deg); Samoa
%     12 = 12 International Date Line (180 deg); Wellington
%     13 = +11 Magadan T. (195 deg); Magadan
%     14 = +10 Australian East S.T. (210 deg); Sydney
%     15 = +9 Japan/Korea S.T (225 deg);Tokyo
%     16 = +8 China S.T. (240 deg); Hong Kong
%     17 = +7 Indochina T. (255 deg); Bangkok
%     18 = +6 Bangladesh T. (270 deg); Omsk
%     19 = +5 Pakistan T. (285 deg); Karachi
%     20 = +4 Samara T. (300 deg); Dubai
%     21 = 0 = +3 Moscow T. (315 deg); Moskau
%     22 = +2 East European T. (330 deg); Kairo
%     23 = +1 Central European T. (345 deg); Berlin
%
% Outputs:
% 
%     Results.TotalRadiation -> total radiation (W/m)
%     Results.DiffuseRadation -> diffuse radiation (W/m)
%     Results.DirectRadiation -> direct radiation (W/m) (multipled by cosine of angle) 
%     Results.ReflectedRadiation -> reflected radiation (W/m)
%     Results.NormalDirectRadiation -> normal direct radiation (W/m)
%     Results.SunAzimuthAngle -> sun's azimuth angle (deg)(tracked power output) 
%     Results.SunElevationAngle -> sun's elevation [altitude] angle (deg)
%     Results.SunIncidenceAngle -> sun's incidence angle (deg)
%     Results.Time -> time vector (minutes)
%
% Example:
% 
%     Time.Hours = 0;
%     Time.Minutes = 0;
%     Time.DayOfYear = 1;
%     Time.StandardMeridian = 8;
%     Time.RunLength = 1;
%     Location.Latitude = 44.5647;
%     Location.Longitude = -123.2608;
%     Location.TiltAngle = 0;
%     Location.GroundReflectance = 0.2;
%     Location.ClearnessFactor = 1;
% 
%     Results = SolarRadiation(Time, Location);
%     plot(Results.Time, Results.DirectRadiation)  
%
%
% Written by: Dr. Hai-Yue Han, Inspired Light
%
% Code Originally Created By: 
%   Dr. Alexandre Marconi de Souza da Costa, State University of Maringa,
%   Brazil
%
%   Bastian Keller,  State University of Maringa, Brazil, and 
%   University of Stuttgart.
%
% Code based on paper:
%
% Keller, B. and Costa, A. M. S. (2011), A Matlab GUI for calculating the 
%   solar radiation and shading of surfaces on the earth. Comput. Appl. 
%   Eng. Educ., 19: 161–170. doi: 10.1002/cae.20301


%% Initializes variables
intMinute = 1;
DayStart = Time.DayOfYear;
DayEnd = DayStart + Time.RunLength-1;


%% Creates system coefficients

%Adds 1 extra month for interpolation wrap around
% AExcel = [1202; 1187;1164;1130;1106;1092;1093;1107;1136;1166;1190;1204;1202];
% BExcel = [0.141;0.142;0.149;0.164;0.177;0.185;0.186;0.182;0.165;0.152;0.142;0.141;0.141];
% CExcel = [0.103;0.104;0.109;0.12;0.13;0.137;0.138;0.134;0.121;0.111;0.106;0.103;0.103];

AExcel = [1202; 1202;1202;1202;1202;1202;1202;1202;1202;1202;1202;1202;1202];
BExcel = [0.141;0.141;0.141;0.141;0.141;0.141;0.141;0.141;0.141;0.141;0.141;0.141;0.141];
CExcel = [0.103;0.103;0.103;0.103;0.103;0.103;0.103;0.103;0.103;0.103;0.103;0.103;0.103];


%Interpolates the data between each month
Ainterp = interp(AExcel,30);
Binterp = interp(BExcel,30);
Cinterp = interp(CExcel,30);

%% Generates sun power and angle data
for intDayOfYear = DayStart:DayEnd

    % read the day of the year and calculate N
    X = intDayOfYear;
    Ndeg = (X - 1) * (360 / 365); % N in deg
    N = Ndeg * pi / 180; % N in rad
    % Equation of time
    EOT = 229.2*(0.000075 + 0.001868*cos(N) - 0.032077*sin(N) - 0.014615*cos(2*N) - 0.04089*sin(2*N));
    %set (handles.TEOT, 'String', num2str (EOT));
    % Sun's declination
    delta = 0.3963723 - 22.9132745*cos(N) + 4.0254304*sin(N) - 0.3872050*cos(2*N) + 0.05196728*sin(2*N) - 0.1545267*cos(3*N) + 0.08479777*sin(3*N);
    %set (handles.Tdelta, 'String', num2str (delta));

    loop = 1440; % Variable of the "for-loop". The loop repeates only one time
    time = 1;
    % Daylight savings time
    hours = Time.Hours;
    minutes = Time.Minutes;
    DST = hours * 60 + minutes; % in minutes


    % Read in location information
    Ll = Location.Longitude; % Longitude
    Ls = (Time.StandardMeridian - 1) * 15; % Standart Meridian
    l = Location.Latitude;   % Latitude

    % Local solar time in minutes
    Lokal_standard_time = DST;
    LST = Lokal_standard_time - (Ll - Ls) * 4 + EOT;
        % Calculation of the hours and minuted for printing out
    if LST < 0
        LSTnegative = 1440 + LST;
        LSTh = floor(LSTnegative/60);
        LSTm = round(rem(LSTnegative,60));
    else
        LSTh = floor(LST/60);
        LSTm = round(rem(LST,60));
    end
    if LSTm == 60
        LSTh = LSTh + 1;
        LSTm = 00;
    end
    % Read surface information
    psi = 0; %str2num (get (handles.EPsi, 'String'));
    alfa = Location.TiltAngle;

    roh = Location.GroundReflectance;
    Cn = Location.ClearnessFactor;

    A = Ainterp(mod(DayStart+intDayOfYear,365)+1);
    B = Binterp(mod(DayStart+intDayOfYear,365)+1);
    C = Cinterp(mod(DayStart+intDayOfYear,365)+1);
    
    
    % Calculation Radiation for 24 hours (loop = 1440) or for setted time (loop
    % = 1)
    %if (ata==0 && handles.hata==0)
    for i = 1:loop
        % Hour angle
        %h = (LST - 720) / 4;
        h = (720 - LST) / 4; %h = LST - 12:00 * (15°/60min)
        % Solar altitude angle
        sinbeta = cos(l*pi/180) * cos(h*pi/180) * cos(delta*pi/180) + sin(l*pi/180) * sin(delta*pi/180);
        beta = asin(sinbeta);
        % Solar azimuth angle
        y_asim = sin(h*pi/180)*cos(delta*pi/180);
        %y_asimtest(i) = y_asim*180;
        z_asim = sin(delta*pi/180)*cos(l*pi/180) - sin(l*pi/180)*cos(h*pi/180)*cos(delta*pi/180);
        %z_asimtest(i) = z_asim*180;
        tanphi = y_asim/z_asim;
        phi = atan(tanphi);
            % check the quadrant of phi
        if z_asim < 0
            phi = phi + pi;
        end
        if phi < 0
            phi = phi + 2*pi;
        end
        % surface solar azimuth
        gama = phi - psi*pi/180;
        % Angle of incidence
        costeta = cos(beta) * cos(gama) * sin(alfa*pi/180) + sinbeta * cos(alfa*pi/180);
        teta = acos(costeta);
        % Calculation of the radiation
        if beta > 0 && beta < pi
            % Normal direct radiation
            Gnd = A*Cn/(exp(B/sinbeta));
            % Konfiguration or angle factor from wall to ground
            Fwg = (1-cos(alfa*pi/180))/2;
            % Direct Radiation
            Gdir = max(costeta, 0)*Gnd;
            % Reflected Radiation
            Gref = roh*Fwg*(sinbeta + C)*Gnd;
            % Diffuse Radiation for Horisontal Surfaces
            GdifH = C*Gnd;
            % Total Radiation (Diffuse Radiation) for tilted surfaces
                % For a nonvertical surface use Fws and for a vertical Gdv/Gdh
            if alfa > 0 && alfa < 90
                Fws = (1+cos(alfa*pi/180))/2;
                Gdif = Fws*GdifH;
                Gtot = Gdir + Gdif + Gref;
                Fss = Fws;
                % Total Radiation (Diffuse Radiation) for vertical surfaces
            elseif alfa == 90
                if costeta > -0.2
                    Gdv_h = 0.55 + 0.437*costeta + 0.313*costeta*costeta;
                else
                    Gdv_h = 0.45;
                end
                Gdif = Gdv_h*GdifH;
                Gtot = Gdir + Gdif + Gref;
                Fss = Gdv_h;
            % Total Radiation for horizontal surfaces
            elseif alfa == 0
                Gdif = GdifH;
                Gtot = Gdir + Gdif; 
                Fss = 1;
            end
        else
            Fss = 0;
            Gnd = 0;
            GdifH = 0;
            Gdir = 0;
            Gdif = 0;
            Gref = 0;
            Gtot = 0;
        end
        % Rise LST one minute, for calculation of 24 hours (loop=1440)
        LST = LST + 1;

        % fill result-vectors 
        Results.TotalRadiation((intDayOfYear-1)*1440+i) = Gtot;
        Results.DiffuseRadation((intDayOfYear-1)*1440+i) = Gdif;
        Results.DirectRadiation((intDayOfYear-1)*1440+i) = Gdir;
        Results.ReflectedRadiation((intDayOfYear-1)*1440+i) = Gref;
        Results.NormalDirectRadiation((intDayOfYear-1)*1440+i) = Gnd;

        Results.SunAzimuthAngle((intDayOfYear-1)*1440+i) = phi*180/pi;
        Results.SunElevationAngle((intDayOfYear-1)*1440+i) = beta*180/pi;
        Results.SunIncidenceAngle((intDayOfYear-1)*1440+i) = teta*180/pi;
        Results.Time((intDayOfYear-1)*1440+i) = intMinute;

        %Increments minutes
        intMinute = intMinute + 1;

        %gamatest(i) = gama*180/pi;
        Fssvect((intDayOfYear-1)*1440+i) = Fss;
        % 3D-PLot
        r = 1;
            % Vectors of the location of the sun in the sky
        threeDx((intDayOfYear-1)*1440+i) = r*cos(beta)*cos(phi);
        threeDy((intDayOfYear-1)*1440+i) = r*cos(beta)*sin(phi);
        threeDz((intDayOfYear-1)*1440+i) = r*sin(beta);
            % Differensation of day and night for the 3D-plot
        if beta > 0
            threeDxa((intDayOfYear-1)*1440+i) = r*cos(beta)*cos(phi);
            threeDya((intDayOfYear-1)*1440+i) = r*cos(beta)*sin(phi);
            threeDza((intDayOfYear-1)*1440+i) = r*sin(beta);
            threeDxb((intDayOfYear-1)*1440+i) = 0;
            threeDyb((intDayOfYear-1)*1440+i) = 0;
            threeDzb((intDayOfYear-1)*1440+i) = 0;
        else
            threeDxb((intDayOfYear-1)*1440+i) = r*cos(beta)*cos(phi);
            threeDyb((intDayOfYear-1)*1440+i) = r*cos(beta)*sin(phi);
            threeDzb((intDayOfYear-1)*1440+i) = r*sin(beta);
            threeDxa((intDayOfYear-1)*1440+i) = 0;
            threeDya((intDayOfYear-1)*1440+i) = 0;
            threeDza((intDayOfYear-1)*1440+i) = 0;
        end

    end
end

end