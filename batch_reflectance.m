%Written by N. Stein 5/30/18
%Perform batch conversion of raw Nano Hyperspec images to reflectance
%INPUTS: raw data cube, header file for raw cube, Nano tarp spectrum, ASD
%tarp spectrum, dark noise spectrum
%NEED TO READ IN AT LEAST TWO FILES FOR THIS CODE TO WORK
%STEPS TO TAKE BEFOREHAND FOR INPUTS. 
%1. Find the tarp in one of the images to get the calibration spectrum. Use
%the brightest segment that isn't over-exposed.
%2. Add .img extension onto data cubes saved out by spectrometer.
%3. Remove all semicolons from the .hdr files because they confuse envihdrread.m. May automate this capability
%later. IF YOU GET "ERROR: INVALID EXPRESSION. CHECK FOR MISSING OR EXTRA CHARACTERS" THIS IS WHY.
%OUTPUTS: reflectance cube, header file for reflectance cube

%READ RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read in airborne spectrometer images
[files,paths] = uigetfile('*.img','MultiSelect','on');

%Read in calibration data
[tarp_asd,asd_txt,asd_raw] = xlsread('C:\Users\ntste\Documents\Ambergris\2017 Data\Spectroscopy\Calibration_Spectra\ambergris_calibration_hyperspec_wavelengths.xlsx'); %Read ASD measurements of tarp spectrum. Parent file doesn't change. These were interpolated to the Nano-Hyperspec wavelengths.
asd_wavelengths = tarp_asd(:,1); %ASD wavelengths
lab_tarp = tarp_asd(:,5); %ASD spectrum from the lightest part of the tarp
field_tarp = tarp_asd(:,4); %ASD spectrum from the medium brightness part of the tarp
%field_tarp = tarp_asd(:,3); %ASD spectrum from the darkest part of the tarp
%cardboard_spectrum = tarp_asd(:,2); %ASD spectrum of the red cardboard
[tarp_raw,tarp_txt,tarp_raw] = xlsread('C:\Users\ntste\Documents\Ambergris\2017 Data\Spectroscopy\Calibration_Spectra\8_3_FLIGHT3\cal_tarp_whitest_5x5.xlsx'); %Read raw tarp spectrum (changes for each scene!)
[dark_raw,dark_txt,dark_raw] = xlsread('C:\Users\ntste\Documents\Ambergris\2017 Data\Spectroscopy\Calibration_Spectra\global_darkref.xlsx'); %Dark spectrum (does not change)
dark_txt(1,:) = []; %Get rid of weird column labels that excel now forces on you
dark_wavelengths = cell2mat(dark_txt(:,1));
dark_wavelengths = str2num(dark_wavelengths); %Convert dark spectrum wavelengths to double
dark_spectrum = cell2mat(dark_txt(:,2));
dark_spectrum = str2num(dark_spectrum); %Convert dark spectrum raw numbers to double
tarp_raw(1,:) = []; 
tarp_wavelengths = cell2mat(tarp_raw(:,1));
tarp_wavelengths = str2num(tarp_wavelengths); %Convert raw tarp wavelengths to double
tarp_field = cell2mat(tarp_raw(:,2));
tarp_field = str2num(tarp_field); %Field measurement of tarp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%BATCH CONVERSION TO REFLECtANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q = 1:length(files) %Iterate over all files
    %Read in header info
    filetemp = cell2mat(files(q));
    file = char(filetemp);
    fileimg = file;
    file(length(file)-2:length(file)) = 'hdr';
    filestr = num2str(file);
    info = envihdrread(strcat(paths,filestr));
    %Read in data cube
    raw_cube = multibandread(strcat(paths,fileimg),[info.lines,info.samples,info.bands],'uint16',0,info.interleave,'ieee-le'); %Read in data cube
    
    reflectance = zeros(info.lines,info.samples,info.bands); %Cube containing reflectance data
    %Convert to reflectance
    for i = 1:info.lines
        i
        for j = 1:info.samples
            for k = 1:info.bands
                reflectance(i,j,k) = (raw_cube(i,j,k)-dark_spectrum(k)).*lab_tarp(k)./(tarp_field(k)-dark_spectrum(k)); %R = (RAW - DARK)*(LAB_TARP_REFL)/(FIELD_TARP_RAW - DARK)
            end
        end
    end
    
    %Write out reflectance file and associated header
    info.data_type = 4;
    file(length(file)-3:length(file)) = [];
    filestr = num2str(file);
    envihdrwrite(info,strcat(paths,filestr,'_reflectance.hdr'));
    fileimg = strcat(paths,filestr,'_reflectance.img');
    reflectance = single(reflectance); %Lower file size
    multibandwrite(reflectance,fileimg,'bil');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

