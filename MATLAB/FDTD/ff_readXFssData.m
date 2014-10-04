function [dataExr, dataExi, dataEyr, dataEyi, dataHxr, dataHxi, dataHyr, dataHyi]  = ff_readXFssData(dataFileExr, dataFileExi, dataFileEyr, dataFileEyi, dataFileHxr, dataFileHxi, dataFileHyr, dataFileHyi, intDataXLength, intDataYLength) 
% Load discrete frequency results from XF into Matlab
%a
% Details in Reference Manual
% Appendix D (Appendix of Calculation Engine Output Files)
% D.2 Multi Point Sensors
%
% Modified by Dr. Hai-Yue Han

%% Get data values

%Opens and parses the Exr files if a Exr file name was given
if size(dataFileExr) > 0   
    
    %Opens the file
    fid = fopen( dataFileExr );
   
    %Reads the entire file
    dataExr = fread( fid, [intDataXLength intDataYLength], 'single' );
    
    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataExr = [];
end

%Opens and parses the Exi files if a Exi file name was given
if size(dataFileExi) > 0   
    
    %Opens the file
    fid = fopen( dataFileExi );
    
    %Reads the entire file
    dataExi = fread( fid, [intDataXLength intDataYLength], 'single' );
    
    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataExi = [];
end

%Opens and parses the Eyr files if a Eyr file name was given
if size(dataFileEyr) > 0   
    
    %Opens the file
    fid = fopen( dataFileEyr );
    
    %Reads the entire file
    dataEyr = fread( fid, [intDataXLength intDataYLength], 'single' );

    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataEyr = [];
end



%Opens and parses the Eyi files if a Eyi file name was given
if size(dataFileEyi) > 0   
    
    %Opens the file
    fid = fopen( dataFileEyi );
    
    %Reads the entire file
    dataEyi = fread( fid, [intDataXLength intDataYLength], 'single' );
    
    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataEyi = [];
end

%Opens and parses the Hxr files if a Hxr file name was given
if size(dataFileHxr) > 0   
    
    %Opens the file
    fid = fopen( dataFileHxr );
   
    %Reads the entire file
    dataHxr = fread( fid, [intDataXLength intDataYLength], 'single' );
    
    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataHxr = [];
end

%Opens and parses the Hxi files if a Hxi file name was given
if size(dataFileHxi) > 0   
    
    %Opens the file
    fid = fopen( dataFileHxi );
    
    %Reads the entire file
    dataHxi = fread( fid, [intDataXLength intDataYLength], 'single' );
    
    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataHxi = [];
end

%Opens and parses the Eyr files if a Eyr file name was given
if size(dataFileHyr) > 0   
    
    %Opens the file
    fid = fopen( dataFileHyr );
    
    %Reads the entire file
    dataHyr = fread( fid, [intDataXLength intDataYLength], 'single' );

    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataHyr = [];
end

%Opens and parses the Eyi files if a Eyi file name was given
if size(dataFileHyi) > 0   
    
    %Opens the file
    fid = fopen( dataFileHyi );
    
    %Reads the entire file
    dataHyi = fread( fid, [intDataXLength intDataYLength], 'single' );
    
    %Closes the file
    fclose(fid);
else 
    %Returns a blank matrix
    dataHyi = [];
end
