function UniversalFileConvert10x
%This code written for the purpose of creating a code capeable of
%efficiently intaking and converting imaris data to the standard matlab
%format OPTIMIZED FOR 10X IMAGES

%Zachary Neronha: Modification History
%Version 1.0 09:23 7 August 2017 Basic Functionality
%Version 1.1 11:56 7 August 2017 Scaling option added
%Version 1.2 14:19 7 August 2017 Bug fix for birth/death mode
%Version 2.0 09:33 9 August 2017: Code now optimized for 10x images
%Version 2.1 10:41 9 August 2017: Code can now store an additional matrix
%with all cell positions (including transient cells)
%Version 2.2 11:38 9 August 2017: Now optimized for automated file loading
%Version 3.0 16:26 25 January 2018: Double dataset support

%THIS IS THE MOST RECENT VERSION AS OF 25 JANUARY 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE TO THE USER: While most of the code is designed to run by itself, the
%input file paths present a slightly different issue, and will need to be
%changed manually for each dataset
%Furthermore the user should take care to upload the POSITION and VELOCITY
%files for the 10x experiments and Vesicles_Nuclei_Position and
%Vesicles_Nuclei_Velocity for the 20x experiments

%Basic Variables to Save: X,Y positions and Velocities
%[2:6 19:23 26:30 43:47 50:54 67:71 74:78 91:95]

clearvars
close all
clc

%select an inputfolder
inputfolder = uigetdir('Z:\ENG_BBCancer_Shared','From where would you like to read the data');
%select an output folder
outputfolder = uigetdir('Z:\ENG_BBCancer_Shared','Where would you like to store the converted data?');
%select wells to process
prompt = 'What wells would you like to process?\n';
wells = input(prompt);
%name the outputfiles
prompt = 'Enter generic file name? (w(number).mat will be added automatically)\n (enter inside single quotes)\n';
outputname = input(prompt);
%look into patch data set
prompt = 'Would you like to store an additional matrix including positions of extremely transient cells? (1=yes)\n';
patchind = input(prompt);

% %get dimensions to load
% prompt = 'Enter sequentially in a matrix the column numbers of x position, y position,positionframe trackID position\n,xVelocity y velocity,velocity frame trackID velocity\n';
% dimenstioncheck = input(prompt);
dimensioncheck = [1 2 6 7 1 2 6 7];
%Xposition Yposition framenumber(pos) trackID(pos) Xvel Yvel
%framenumber(vel) trackID(vel)

%If rescaling is needed determine the ranges
prompt = 'Would you like to rescale the position data? (1(yes) 0 (no))\n';
rescaleID = input(prompt);
if rescaleID == 1
%    prompt = 'Input final x width\n';
%    xbounds = input(prompt);
%    prompt = 'Input final y height\n';
%    ybounds = input(prompt);    
    xbounds = [382 1282];
    ybounds = [252 1152];
    disp('HARD CODED VALUES USED TO SCALE TO 900x900');
end
tic
for well = wells
    
    %CHANGE THE FILE PATHS HERE IF NECESSARY!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%
    filetoload = strcat(inputfolder,'\w',num2str(well),'_Statistics\w',...
        num2str(well),'_Position.csv');
    posdata = xlsread(filetoload);

    %load the velocity data
    filetoload = strcat(inputfolder,'\w',num2str(well),'_Statistics\w',...
        num2str(well),'_Velocity.csv');
    veldata = xlsread(filetoload);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %account for the birth/death issue
    if posdata(1,7) == 900
        posframes = round(posdata(:,dimensioncheck(3)+1)./900);
        velframes = round(veldata(:,dimensioncheck(7)+1)./900);
        Xposition = posdata(:,dimensioncheck(1));
        Yposition = posdata(:,dimensioncheck(2));        
        posTrackID = posdata(:,dimensioncheck(4)+1);
        Xvel = veldata(:,dimensioncheck(5));
        Yvel = veldata(:,dimensioncheck(6));
        velTrackID = veldata(:,dimensioncheck(8)+1);
        disp('Converted from Birth Death Mode');
    else
        %relabel and extract as directed by the user
        Xposition = posdata(:,dimensioncheck(1));
        Yposition = posdata(:,dimensioncheck(2));
        posframes = posdata(:,dimensioncheck(3));
        posTrackID = posdata(:,dimensioncheck(4));
        Xvel = veldata(:,dimensioncheck(5));
        Yvel = veldata(:,dimensioncheck(6));
        velframes = veldata(:,dimensioncheck(7));
        velTrackID = veldata(:,dimensioncheck(8));

    end

    %rescale the data
    posTrackID = (posTrackID-min(posTrackID))+1;
    velTrackID = (velTrackID-min(velTrackID))+1;
    Xposition = Xposition - min(Xposition);
    Yposition = Yposition - min(Yposition);
    
    if patchind == 1
        %store all positions
        patchX = Xposition;
        patchY = Yposition;
        savename2 = strcat(outputfolder,'\',outputname,'AllCells_w',num2str(well),'.mat');
        if rescaleID == 1
            StoreAllPositions(patchX,patchY,xbounds,ybounds,rescaleID,posTrackID,posframes,savename2);
        else
            StoreAllPositions(patchX,patchY,0,0,rescaleID,posTrackID,posframes,savename2);
        end
    end
    
    %get rid of extremely transient cells
    m = isnan(posTrackID) == 0; 
    Xposition = Xposition(m);
    Yposition = Yposition(m);
    posframes = posframes(m);
    posTrackID = posTrackID(m);
    
    if rescaleID == 1
        %determine if the x value violates our bounds
        m = (Xposition < xbounds(1))|(Xposition > xbounds(2))|...
            (Yposition < ybounds(1))|(Yposition > ybounds(2));
        m = ~m; %convert logical operator
        %keep only the values that don't violate our standards
        Xposition = Xposition(m);
        Yposition = Yposition(m);
        Xvel = Xvel(m);
        Yvel = Yvel(m);
        posTrackID = posTrackID(m);
        velTrackID = velTrackID(m);
        posframes = posframes(m);
        velframes = velframes(m);
        
        %rescale the data one final time
        Xposition = Xposition - xbounds(1);
        Yposition = Yposition - ybounds(1);
    end
    
    %preallocate storage matricies for operation speed
    storeX = NaN(max(posTrackID),max(posframes));
    storeY = NaN(max(posTrackID),max(posframes));
    storevelX = NaN(max(velTrackID),max(velframes));
    storevelY = NaN(max(velTrackID),max(velframes));
    

   %store the data in the appropriate location
    for columnloop = 1:size(Xposition,1)
       storeX(posTrackID(columnloop),posframes(columnloop)) = Xposition(columnloop); 
       storeY(posTrackID(columnloop),posframes(columnloop)) = Yposition(columnloop); 
       storevelX(velTrackID(columnloop),velframes(columnloop)) = Xvel(columnloop); 
       storevelY(velTrackID(columnloop),velframes(columnloop)) = Yvel(columnloop); 
       u = isnan(Xposition(columnloop))+isnan(Yposition(columnloop))+isnan(Xvel(columnloop))+isnan(Yvel(columnloop));
       
    end

    savename = strcat(outputfolder,'\',outputname,'w',num2str(well),'.mat');
    save(savename,'storeX','storeY','storevelX','storevelY');
    fprintf('Well %d is now complete!\n',well);
    beep
end

toc
disp('CONVERSION COMPLETE');
end

function StoreAllPositions(patchX,patchY,xbounds,ybounds,rescaleID,posTrackID,posframes,savename)
    %rescale the data
    if rescaleID == 1
        %determine if the x value violates our bounds
        m = (patchX < xbounds(1))|(patchX > xbounds(2))|...
            (patchY < ybounds(1))|(patchY > ybounds(2));
        m = ~m; %convert logical operator
        %keep only the values that don't violate our standards
        patchX = patchX(m);
        patchY = patchY(m);
        posTrackID = posTrackID(m);
        posframes = posframes(m);

        %rescale the data one final time
        patchX = patchX - xbounds(1);
        patchY = patchY - ybounds(1);
    end
    
    %assign trackIDs to the cells that lack them
    TID = max(posTrackID) + 1;
    for aa = 1:size(patchX,1)
        if isnan(posTrackID(aa)) == 1
            posTrackID(aa) = TID;
            TID = TID + 1;
        end
    end
    %preallocate matricies to store data in
    PstoreX = nan(max(posTrackID),max(posframes));
    PstoreY = nan(max(posTrackID),max(posframes));

    %store the data in the appropriate location
    for columnloop = 1:size(patchX,1)
       PstoreX(posTrackID(columnloop),posframes(columnloop)) = patchX(columnloop); 
       PstoreY(posTrackID(columnloop),posframes(columnloop)) = patchY(columnloop); 
    end
    save(savename,'PstoreX','PstoreY');
end