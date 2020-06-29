%Written by N. Stein 6/11/18
%Reads in IMU and GPS data, converts to UTM, and determines the location of
%the upper left pixel in UTM space
%N up, S down --> no rotation, first frame is top
%N down, S up--> rotate 180 degrees, first frame is now at bottom, so use last frame for location

%Read in imu data and convert to utm
[num1,txt1,raw1] = xlsread('C:\nanoImgs\8_3_FLIGHT3\imu_gps.xlsx');
altitudes = num1(:,6);
zero_altitude = 34.358; %global datum means ground level is negative, needs to be added for relative altitude
cols = 640.0; %number of columns
[x,y,utmzone] = ll2utm([num1(:,4),num1(:,5)],'nad27');

%Find the orientation of the scene
[num,txt,raw] = xlsread('C:\nanoImgs\8_3_FLIGHT3\frameIndex_0.xlsx');
frame_time = num(:,2);
for i = 1:length(frame_time)
    time_diff = abs(num1(:,7)-frame_time(i));
    row_index_check = find(time_diff==min(time_diff));
    row_index(i) = row_index_check(1);
    xloc(i) = x(row_index(i));
    yloc(i) = y(row_index(i));
end

xdiff = xloc(length(frame_time))-xloc(1); %Total change in x location
ydiff = yloc(length(frame_time))-yloc(1); %Total change in y location
fracdiff = xdiff/ydiff; %compare fractional displacement between x and y coordinates
%N-S oriented w/ N up: y has negative slope and fracdiff < 0.05. First frame at top.
%N-S oriented w/ S up: y has positive slope and fracdiff < 0.05. First frame at bottom. Rotate 180 deg.

%Perform scene rotation if necessary
if (ydiff < 0) && (fracdiff < 0.05)
    'N-S w/ N up. No reorientation necessary.'
    a = multibandread('C:\nanoImgs\8_3_FLIGHT3\raw_8080_reflectance_bd680.img',[8064, 640, 1],'float',0,'bsq','ieee-le');
    a = single(a);
    y_comparison_factor = 1;
    check_1 = 1;
    %x_rotation_factor = 1;
elseif (ydiff > 0) && (fracdiff < 0.05)
    'N-S w/ S up. Rotating by 180 degrees.' %Use last frame for y
    y_comparison_factor = -1;
    %x_rotation_factor = -1;
    a = multibandread('C:\nanoImgs\8_3_FLIGHT3\raw_0_reflectance_bd680.img',[8080, 640, 1],'float',0,'bsq','ieee-le');
    a = rot90(rot90(a));
    a = single(a);
    frame_offset = 0;
    check_1 = 0;
    %multibandwrite(a,'C:\nanoImgs\8_3_FLIGHT3\raw_0_reflectance_bd680_rot.img','bil');
else
    'Expected orientation case not met.'
end

%Find x,y point that corresponds to the left pixel in each row
for i = 1:length(frame_time)
    i
    altitude(i) = zero_altitude + altitudes(row_index(i)); %Use altitude of each row
    resolution_x(i) = altitude(i)*.001428; %x Pixel size in meters, .00155
    resolution_y(i) = resolution_x(i); %No correction for smear
    left_x(i) = xloc(i) - resolution_x(i)*(cols./2.0); %UTM x coord of left-most pixel
    left_y(i) = yloc(i); %UTM y corod of left-most pixel
end

%Save out to files
master_hdr = envihdrread('C:\nanoImgs\8_3_FLIGHT3\raw_0_reflectance_bd680_rot.hdr');
path = 'C:\nanoImgs\8_3_FLIGHT3\row_by_row\';
for q = 1:length(frame_time)
    q
    master_hdr.lines = 1; %only saving one row at a time
    master_hdr.bands = 1; 
    mapinfo = master_hdr.map_info;
    mapinfo = strsplit(mapinfo);
    mapinfo(4) = cellstr(strcat(num2str(left_x(q*check_1 + (1-check_1)*(length(frame_time)-q+1))),',')); %Update x and y UTM info in header
    mapinfo(5) = cellstr(strcat(num2str(left_y(q*check_1 + (1-check_1)*(length(frame_time)-q+1))),','));
    mapinfo(6) = cellstr(strcat(num2str(resolution_x(q*check_1 + (1-check_1)*(length(frame_time)-q+1))),',')); %Update x and y pixel size in header
    mapinfo(7) = cellstr(strcat(num2str(resolution_y(q*check_1 + (1-check_1)*(length(frame_time)-q+1))),',')); %Update x and y pixel size in header
    master_hdr.map_info = strjoin(mapinfo);
    envihdrwrite(master_hdr,strcat(path,'raw_0_reflectance_bd680_',num2str(q),'.hdr'));    
    multibandwrite(a(q,:),strcat(path,'raw_0_reflectance_bd680_',num2str(q),'.img'),'bil');
end
