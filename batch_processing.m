%Written by N. Stein 5/30/18
%Perform batch conversion of raw Nano Hyperspec images to reflectance. Also
%calculates and save band depths and centers for select pigments.
%INPUTS: raw data cube, header file for raw cube, Nano tarp spectrum, ASD
%tarp spectrum, dark noise spectrum
%NEED TO READ IN AT LEAST TWO FILES FOR THIS CODE TO WORK
%STEPS TO TAKE BEFOREHAND FOR INPUTS. 
%1. Find the tarp in one of the images to get the calibration spectrum. Use
%the brightest segment that isn't over-exposed.
%2. Add .img extension onto data cubes saved out by spectrometer.
%3. Remove all semicolons from the .hdr files because they confuse envihdrread.m. May automate this capability
%later. IF YOU GET "ERROR: INVALID EXPRESSION. CHECK FOR MISSING OR EXTRA CHARACTERS" THIS IS WHY.
%OUTPUTS: reflectance cube, header file for reflectance cube, band depth
%and band center maps.

%READ RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read in airborne spectrometer images
[files,paths] = uigetfile('*.img','MultiSelect','on');

%Read in calibration data
[tarp_asd,asd_txt,asd_raw] = xlsread('C:\Users\ntste\Documents\Ambergris\2017 Data\Spectroscopy\Calibration_Spectra\ambergris_calibration_hyperspec_wavelengths.xlsx'); %Read ASD measurements of tarp spectrum. Parent file doesn't change. These were interpolated to the Nano-Hyperspec wavelengths.
asd_wavelengths = tarp_asd(:,1); %ASD wavelengths
lab_tarp = tarp_asd(:,5); %ASD spectrum from the lightest part of the tarp
%field_tarp = tarp_asd(:,4); %ASD spectrum from the medium brightness part of the tarp
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
%tarp_wavelengths = cell2mat(tarp_raw(:,1));
%tarp_wavelengths = str2num(tarp_wavelengths); %Convert raw tarp wavelengths to double
tarp_field = cell2mat(tarp_raw(:,2));
%tarp_field = str2num(tarp_field); %Field measurement of tarp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SET CONSTANTS FOR PARAMETER MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pigments_to_check = ['590','630','680','720','740','790','800','840','870','900']; %String of bands to check (nm)
%Chl a
start_index_590 = 82;
end_index_590 = 106;
acceptable_left_590 = 83;
acceptable_right_590 = 91;
start_index_630 = 83;
end_index_630 = 115;
acceptable_left_630 = 87;
acceptable_right_630 = 108;
start_index_680 = 114;
end_index_680 = 141;
acceptable_left_680 = 124;
acceptable_right_680 = 132;
start_index_720 = 143;
end_index_720 = 151;
acceptable_left_720 = 144;
acceptable_right_720 = 150;
start_index_740 = 133;
end_index_740 = 170;
acceptable_left_740 = 150;
acceptable_right_740 = 161;
start_index_790 = 150;
end_index_790 = 189;
acceptable_left_790 = 173;
acceptable_right_790 = 179;
start_index_800 = 164;
end_index_800 = 190;
acceptable_left_800 = 179;
acceptable_right_800 = 182;
start_index_840 = 190;
end_index_840 = 219;
acceptable_left_840 = 197;
acceptable_right_840 = 206;
start_index_870 = 193;
end_index_870 = 238;
acceptable_left_870 = 208;
acceptable_right_870 = 216;
start_index_900 = 217;
end_index_900 = 249;
acceptable_left_900 = 222;
acceptable_right_900 = 231;
%
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
    
    %MAKE PARAMETER MAPS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate 680 nm BD and BC (chl. a)
    counter = 1;
    while counter < length(pigments_to_check)
        %Assign indices to the relevant bands
        band_string = pigments_to_check(counter:counter+2);
        start_index = eval(strcat('start_index_',band_string));
        end_index = eval(strcat('end_index_',band_string));
        acceptable_left = eval(strcat('acceptable_left_',band_string));
        acceptable_right = eval(strcat('acceptable_right_',band_string));
        sample_res =start_index:end_index;
        wavelengths = asd_wavelengths(sample_res);
        accept_res = acceptable_left:acceptable_right;
        for l = 1:info.lines
            disp(['line ',num2str(l),' of ',num2str(info.lines),' in BD ',num2str(band_string)])
            for m = 1:info.samples
                y = squeeze(reflectance(l,m,:));
                y = smooth(y); %denoise with 5 point running average
                loc_maxpoint = start_index + find(y(start_index+1:end_index-1) == max(y(start_index+1:end_index-1))); %Locate the maximum point between the initial indices
                nodes = [start_index,loc_maxpoint(1),end_index];
                nodecheck = nodes(1:3);
                nodecheck = unique(nodecheck); %No duplicate nodes
                complete = 0;

                while complete == 0
                    newnode_counter = 0;
                    newnodes = [];
                    vals = [];
                    for i = 1:length(nodecheck)-1
                        slope = (y(nodecheck(i+1)) - y(nodecheck(i)))/(nodecheck(i+1) - nodecheck(i));          %Slope between nodes
                        line = 0:(nodecheck(i+1)-nodecheck(i));
                        line = y(nodecheck(i)) + line*slope;
                        vals = [y(nodecheck(i):nodecheck(i+1))];  
                        div = vals./line';                                                                      %Divide reflectance spectrum by line between nodes
                        divgt = div(div > 1.001);
                        if sum(divgt) > 0 && (nodecheck(i+1) - nodecheck(i)) > 1                              %If div > 1 anywhere along segment then add new node. Do not investigate adjacent nodes.
                            newnode_counter = newnode_counter + 1;
                            nodepos = find(div == max(divgt));
                            newnodes(newnode_counter) = -1 + nodecheck(i) + nodepos(1); %Calculate position of new node as max in segment          
                            if newnodes(1) == end_index
                                newnode_counter = [];
                                newnodes = [];
                            end
                        end
                    end
                    complete = 1 - any(newnode_counter);                              %End if no new nodes are added
                    nodecheck = [nodecheck newnodes];                                 %Add new nodes
                    nodecheck = sort(nodecheck);                                      %Sort nodes
                    nodecheck = unique(nodecheck);                                    %No duplicate nodes
                end
                hull_values = interp1(nodecheck,y(nodecheck),sample_res);
                cont_rem_spec = y(nodecheck(1):nodecheck(length(nodecheck)))-hull_values'; %Continuum removed spectrum
                cont_rem_spec_subset = cont_rem_spec(acceptable_left-start_index:acceptable_right-start_index);
                v1 = -1.0*min(cont_rem_spec);
                v2 = wavelengths(find(cont_rem_spec == min(cont_rem_spec)));
                bd(l,m) = v1(1);                                      %Band depth (positive value)
                bc(l,m) = v2(1);       %Band center
            end
        end
        multibandwrite(single(bd),strcat(paths,filestr,'_reflectance_bd',band_string,'.img'),'bil');
        multibandwrite(single(bc),strcat(paths,filestr,'_reflectance_bc',band_string,'.img'),'bil');
        info.bands = 1;
        envihdrwrite(info,strcat(paths,filestr,'_reflectance_bd',band_string,'.hdr'));
        envihdrwrite(info,strcat(paths,filestr,'_reflectance_bc',band_string,'.hdr'));
        counter = counter + 3; %Iterating by 3 because each band wavelength is a three digit number
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

