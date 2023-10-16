%% Initialise
% Extensometer with MTS data and mp4 files 
% Version15

% Clear figures and variables
clear
close all

% Select video folder
Names=dir('data\vids');
Slash='\';

% 1= Show progression of ROI during strain measurement
% 0= Do not show progression (takes less time)
Show=0;

% Machine Sample rate
MTSHz=102.4;

% CameraFPS
FPS=60;

% Set Level of Otsu Filter
SetLevel=0.3;

% Strain Points

StrainStart=0.002;
StrainEnd=0.004;

% Insert diameters in order of sample no. for stress calculations.

Diameters=[5]; 

% Number of points per segment of best fit line

SegmentLength=150;

% Polynomial degree for best fit line (1 or 2)

Degree=1;

%Stress Threshold - below these values will not be used for Young's Modulus
%calculations

StressThreshold=0.15; %MPa

% Line size is the size of black row needed to constitute an edge

LineSize=50;

% Number of white pixels that does not interrupt the location of edge top
% and bottom

WhiteLine=10;

% Move points at which modulus is calculated if strain does not start at 0 (0=no, 1=yes)
AdjustModulusPoints=0;

% If want a single curve that is nicely smoothed but less accurate, 1, if not 0.
SmoothInaccurate=0;

% Plot Raw Data 

RawData=0;

%% Get Frames

% Cut out junk directories

Counter2=1;
for i=1:length(Names)
    Name=(Names(i,1));
    Temp=double(Name.isdir);
    if Temp==1
        Counter2=Counter2+1;
    end
end

% Put Names of videos into an array
for i= Counter2:length(Names)
    Name=(Names(i,1));
    NameTxt=(Name.name);
    NameArray{i}=NameTxt;
end 

NameArray=NameArray(Counter2:end);

FolderTxt=Names.folder;

% Create DataCell to store all data

DataCell=cell((length(NameArray)),4);


for i=(1:length(NameArray))
    
% Choose videos one at a time
FullLocation=strcat(FolderTxt,Slash,NameArray{i});  
    
% Open video and find number of frames
VideoFile=VideoReader(FullLocation);
FramesNo=VideoFile.Duration*VideoFile.FrameRate;
FramesNo=floor(FramesNo);
FrameArray={};

% Insert Name into DataCell
DataCell{i,1}=NameArray{i};

counter=0;
% Remove all frames from video and ensure they are all numbered
while hasFrame(VideoFile)
counter=counter+1;
number=num2str(counter);
buffer='00000';

if (counter >= 10)
    buffer='0000';
end
if (counter<=100) && (counter >= 10)
    buffer='000';
end
if (counter<=1000) && (counter>=100) && (counter >= 10)
    buffer='00';
end
if (counter<=10000) && (counter>=1000) && (counter>=100) && (counter >= 10)
    buffer='0';
end
% Convert frames to Otsu images
IndividualFrame=readFrame(VideoFile);
IndividualFrame=rgb2gray(IndividualFrame);
% Level used in Otsu filter can be changed
level=SetLevel;
IndividualFrame=imbinarize(IndividualFrame,level); 
% Rotate frame for processing
IndividualFrame=imrotate(IndividualFrame,90);
% Create cell with all frames of video
if counter==1
    FrameArray={IndividualFrame};
else
    FrameArray{end+1}=IndividualFrame;
end

end
% Add frame cell to main DataCell
  DataCell{i,2}=FrameArray;
end

%% Initial Conditions

% To increase reproducibility of results, automate the ROI selection

for IMAGE=(1:length(NameArray))
 FrameLength=length(DataCell{IMAGE,2});
 
% INPUT MIDDLE, TOP, BOTTOM (APPROX)
FirstImStart=DataCell{IMAGE,2};
FirstIm=FirstImStart(1);
FirstIm=cell2mat(FirstIm);
imshow(FirstIm)
[x1,y1]=ginput(3);

%Take middle X coordinate, as well as top and bottom y positions

MiddleROI=floor(x1(1,1));
BottomROI1=floor(y1(3,1));
TopROI1=floor(y1(2,1));

%LEFT SIDE

% Flags used later
flag1=0;
flag4=0;
flag5=0;
%Count from middle towards left
for i=(MiddleROI:-1:1)
    % Select each row moving left via "i"
    ValueArray=FirstIm(TopROI1:BottomROI1,i);
    Counter1_initial=0;
    
    % Run through the values in each row
    for j=(1:length(ValueArray))
        % Only run if ROI's have not been established
        if flag4==0
            Value=ValueArray(j);
            % If pixel is black and line has not been established, add to
            % count
            if Value==0 && Counter1_initial<LineSize && flag1==0
                Counter1_initial=Counter1_initial+1;
            end
            % If pixel if white and line has not been established, reset
            if Value==1 && Counter1_initial<LineSize && flag1==0
                Counter1_initial=0;
            end
            % If line has been established and pixel is white, set leftROI,
            if Counter1_initial >= LineSize && flag1==0 
                LeftROI=i;
                TestROI=LeftROI;
                % All if loops nullified
                flag1=1;
            end
        end
    end
    % If LeftROI has been founded and this loop has not yet run, run it.
    if flag1==1 && flag4==0
        % Initialise counters
        flag3=0;
        Black_initial=0;
        White_initial=0;
        
        % Select LeftROI row
        Row=FirstIm(TopROI1:BottomROI1,LeftROI);
        
        % Run through each value of the row
        for j=(1:length(Row))
            Value=Row(j);
            % For first black pixel, set bottomROI
            if Value==0 && flag3==0 && Black_initial<LineSize
                flag3=1;
                TopROI_Left=j+TopROI1;
                Black_initial=Black_initial+1;
            end
            
            % For subsequent black pixels, add one to counter
            if Value==0 && flag3==1 && Black_initial<LineSize
                Black_initial=Black_initial+1;
                % If white pixels have been counted, reset this counter
                White_initial=0;
            end
            
            % If black line has not been established and white pixel
            % appears, bin the line and start again
            if Value==1 && flag3==1 && Black_initial<LineSize && flag5==0
                flag3=0;
                Black_initial=0;
            end
            
            % If pixel is white and black line has been established, start
            % white counter
            if Value==1 && flag3==1 && Black_initial>=LineSize && White_initial<WhiteLine
                White_initial=White_initial+1;
                Black_initial=Black_initial+1;
            end
            
            % If pixel is still white after 10 counts, decide the line
            % ended at the first white pixel registered and save topROI
            if Value==1 && flag3==1 && Black_initial>=LineSize &&  White_initial>=WhiteLine
                Black_initial=Black_initial-White_initial;
                BottomROI_Left=TopROI_Left+Black_initial;
                flag5=1;
            end
        end
        % Don't run once vertical ROI's are found
        flag4=1;  
    end
end

%RIGHT SIDE

% Flags used later
flag1=0;
flag4=0;
flag5=0;
%Count from middle towards rigth
for i=(MiddleROI:length(FirstIm))
    % Select each row moving left via "i"
    ValueArray=FirstIm(TopROI1:BottomROI1,i);
    Counter1_initial=0;
    
    % Run through the values in each row
    for j=(1:length(ValueArray))
        % Only run if ROI's have not been established
        if flag4==0
            Value=ValueArray(j);
            % If pixel is black and line has not been established, add to
            % count
            if Value==0 && Counter1_initial<LineSize && flag1==0
                Counter1_initial=Counter1_initial+1;
            end
            % If pixel if white and line has not been established, reset
            if Value==1 && Counter1_initial<LineSize && flag1==0
                Counter1_initial=0;
            end
            % If line has been established and pixel is white, set letROU,
            if Counter1_initial >= LineSize && flag1==0 
                RightROI=i;
                % All if loops nullified
                flag1=1;
            end
        end
    end
    % If RightROI has been founded and this loop has not yet run, run it.
    if flag1==1 && flag4==0
        % Initialise counters
        flag3=0;
        Black_initial=0;
        White_initial=0;
        
        % Select LeftROI row
        Row=FirstIm(TopROI1:BottomROI1,RightROI);
        
        % Run through each value of the row
        for j=(1:length(Row))
            Value=Row(j);
            
            % For first black pixel, set bottomROI
            if Value==0 && flag3==0 && Black_initial<LineSize
                flag3=1;
                TopROI_Right=j+TopROI1;
                Black_initial=Black_initial+1;
            end
            
            % For subsequent black pixels, add one to counter
            if Value==0 && flag3==1 && Black_initial<LineSize
                Black_initial=Black_initial+1;
                % If white pixels have been counted, reset this counter
                White_initial=0;
            end
            
            % If black line has not been established and white pixel
            % appears, bin the line and start again
            if Value==1 && flag3==1 && Black_initial<LineSize && flag5==0
                flag3=0;
                Black_initial=0;
            end
            
            % If pixel is white and black line has been established, start
            % white counter
            if Value==1 && flag3==1 && Black_initial>=LineSize && White_initial<WhiteLine
                White_initial=White_initial+1;
                Black_initial=Black_initial+1;
            end
            
            % If pixel is still white after 10 counts, decide the line
            % ended at the first white pixel registered and save topROI
            if Value==1 && flag3==1 && Black_initial>=LineSize &&  White_initial>=WhiteLine
                Black_initial=Black_initial-White_initial;
                BottomROI_Right=TopROI_Right+Black_initial;
                flag5=1;
            end
        end
        % Don't run once vertical ROI's are found
        flag4=1;  
    end
end

% Decide which side has a higher TopROI and which side  has a lower BottomROI
TempTop=TopROI_Left-TopROI_Right;
TempBottom=BottomROI_Left-BottomROI_Right;
if TempTop<1
    TopROI=TopROI_Right;
else
    TopROI=TopROI_Left;
end

if TempBottom<1
   BottomROI=BottomROI_Left;
else
    BottomROI=BottomROI_Right;
end

StrainArray=zeros(1,FrameLength);

% Cut Image
FirstImCut=FirstIm(TopROI:BottomROI,LeftROI:RightROI);
imshow(FirstImCut)

MidPoint=(length(FirstImCut)/2);
MidFull=floor(LeftROI+MidPoint);

FirstImCut=FirstIm(TopROI:BottomROI,:);

%Define Initial length
InitialLength=RightROI-LeftROI;


%% Record expansion

% Set Initial Variables
Expansion=0;

LeftEdge=zeros(1,FrameLength);
RightEdge=zeros(1,FrameLength);
ImageStrain=0;
for i=(1:FrameLength)

    CounterLeft=0;
    CounterRight=0;
    flagright=0;
    flagleft=0;
    if i==1
        for j=(MidFull:-1:1)
            Row=FirstImCut(:,j);
            [vert,horz]=size(FirstImCut);
            for k=(1:length(Row))
                Value=Row(k);
                if flagleft==0
                    if Value == 1 && CounterLeft<(vert)
                        CounterLeft=0;
                    end
                    if Value==0 && CounterLeft<(vert-1)
                        CounterLeft=CounterLeft+1;
                    end
                    if CounterLeft>=(vert-1) 
                        LeftEdge(i)=j;
                        flagleft=1;
                    end
                end
            end
        end
    
        for j=(MidFull:length(FirstImCut))
            Row=FirstImCut(:,j);
            [vert,horz]=size(FirstImCut);
            for k=(1:length(Row))
                Value=Row(k);
                if flagright==0
                    if Value == 1 && CounterRight<(vert)
                        CounterRight=0;
                    end
                    if Value==0 && CounterRight<(vert-1)
                        CounterRight=CounterRight+1;
                    end
                    if CounterRight>=(vert-1) 
                        RightEdge(i)=j;
                        flagright=1;
                    end
                end
            end
        end
        ImageStrain=0;
    else
        
    CurrIm=DataCell{IMAGE,2};
    CurrIm=CurrIm(i);
    CurrIm=cell2mat(CurrIm);
    CurrImCut=CurrIm(TopROI:BottomROI,:);
    leftflag=0;
    rightflag=0;
    
    for j=(MidFull:-1:1)
        Row=CurrImCut(:,j);
        [vert,horz]=size(FirstImCut);
        for k=(1:length(Row))
            Value=Row(k);
            if flagleft==0
                if Value == 1 && CounterLeft<(vert) 
                    CounterLeft=0;
                end
                if Value==0 && CounterLeft<(vert-1)
                    CounterLeft=CounterLeft+1;
                end
                if CounterLeft>=(vert-1) 
                    LeftEdge(i)=j;
                    flagleft=1;
                    LeftDisp=LeftEdge(i)-LeftEdge(i-1);
                end
            end
        end
    end 
        for j=(MidFull:length(FirstImCut))
            Row=FirstImCut(:,j);
            [vert,horz]=size(FirstImCut);
            for k=(1:length(Row))
                Value=Row(k);
                if flagright==0
                    if Value == 1 && CounterRight<(vert)
                        CounterRight=0;
                    end
                    if Value==0 && CounterRight<(vert-1)
                        CounterRight=CounterRight+1;
                    end
                    if CounterRight>=(vert-1) 
                        RightEdge(i)=j;
                        flagright=1;
                        RightDisp=RightEdge(i)-RightEdge(i-1);
                    end
                end
            end
        end
        ImageStrainTemp=-LeftDisp+RightDisp;
        ImageStrain=ImageStrain+ImageStrainTemp;
    end
    % Find Strain and store in array
    StrainArray(i)=ImageStrain/InitialLength;
    if Show==1 && i>1
        imshow(CurrImCut)
    end
end
% Store strain data in DataCell
DataCell{IMAGE,3}=StrainArray;
% Create storage cells

end

StrainCell=cell(1,length(NameArray));

for i=(1:length(NameArray))
    % Update figure counter (ensures nothing plots over each other)
    if i==1
        FigureCounter=1;
    else
        FigureCounter=FigureCounter+1;
    end

%Store number of frames    
FrameLength=length(DataCell{i,2});
%Define number of images in array
ImageNo=linspace(1,FrameLength,FrameLength);
% Get Strain of Image
StrainArray=DataCell{i,3};

%
% Plot strain and update counter
figure(FigureCounter)
plot(ImageNo,StrainArray,'o')
FigureCounter=FigureCounter+1;


BestFitArray=zeros(1,length(ImageNo));

if SmoothInaccurate==0 && RawData==0
% Best fit each segment as defined by segment length and polynomial degree

LoopEnd=ceil(length(StrainArray)/SegmentLength);
for j=(1:LoopEnd)
    if j==1
        StrainArrayCurr=StrainArray(1:j*SegmentLength);
        ImageNoCurr=ImageNo(1:j*SegmentLength);
        BestFit2=polyfit(ImageNoCurr,StrainArrayCurr,Degree);
        for k=(1:j*SegmentLength)
            if Degree==1
            BestFitArray(k)=BestFit2(1)*k+BestFit2(2);
            end
            if Degree==2
            BestFitArray(k)=BestFit2(1)*k^2+BestFit2(2)*k+BestFit2(3);
            end
        end
    end
    if j>1 && j<LoopEnd
        StrainArrayCurr=StrainArray(((j-1)*SegmentLength+1):j*SegmentLength);
        ImageNoCurr=ImageNo(((j-1)*SegmentLength+1):j*SegmentLength);
        BestFit2=polyfit(ImageNoCurr,StrainArrayCurr,Degree);
        for k=(((j-1)*SegmentLength+1):j*SegmentLength)
            if Degree==1
            BestFitArray(k)=BestFit2(1)*k+BestFit2(2);
            end
            if Degree==2
            BestFitArray(k)=BestFit2(1)*k^2+BestFit2(2)*k+BestFit2(3);
            end
        end
    end 
    if j==LoopEnd
        StrainArrayCurr=StrainArray(((j-1)*SegmentLength):end);
        ImageNoCurr=ImageNo(((j-1)*SegmentLength):end);
        BestFit2=polyfit(ImageNoCurr,StrainArrayCurr,Degree);
        for k=(((j-1)*SegmentLength+1):length(ImageNo))
            if Degree==1
            BestFitArray(k)=BestFit2(1)*k+BestFit2(2);
            end
            if Degree==2
            BestFitArray(k)=BestFit2(1)*k^2+BestFit2(2)*k+BestFit2(3);
            end
        end
    end
end
end

if SmoothInaccurate==1 && RawData==0
    BestFitSmooth=polyfit(ImageNo,StrainArray,2);
    for j=(1:length(StrainArray))
        BestFitArray(j)=BestFitSmooth(1)*j^2+BestFitSmooth(2)*j+BestFitSmooth(3);
    end
end

if RawData==1
    BestFitArray=StrainArray;
end
%Plot Strain with best fit line
figure(FigureCounter)

plot(ImageNo,StrainArray,'o')
hold on
plot(ImageNo,BestFitArray)
hold off
% Store Best fit
StrainCell{1,i}=BestFitArray;
end

%% Get Stress

% Take stress data from folder
Names=dir('data\TensileData');
NameArrayStress=cell((length(Names)),1);

Counter2=1;
for i=1:length(Names)
    Name=(Names(i,1));
    Temp=double(Name.isdir);
    if Temp==1
        Counter2=Counter2+1;
    end
end

for i= Counter2:length(Names)
    Name=(Names(i,1));
    NameTxt=(Name.name);
    NameCells{i}=NameTxt;
end 

Counter2=1;
for i=1:length(Names)
    Name=(Names(i,1));
    Temp=double(Name.isdir);
    if Temp==1
        Counter2=Counter2+1;
    end
end

% Insert Names into Array
for i= Counter2:length(Names)
    Name=(Names(i,1));
    NameTxt=(Name.name);
    NameArrayStress{i}=NameTxt;
end 

%Cut out junk
NameArrayStress=NameArrayStress(Counter2:end);

for i=(1:length(NameArrayStress))
%Create folder directory for locating data
Folder=(Names(i,1));
FolderTxt=(Folder.folder);
Slash='\';
FolderTxt=strcat(FolderTxt,Slash);

Location=strcat(FolderTxt,NameArrayStress{i});
    ID=fopen(Location,'r');
    % %s's are for each column of data and they are imported as strings
    Data=textscan(ID,'%s %s %s %s %s');
    fclose(ID);
    
    ArrayDimensions=zeros(length(Data{1,1}),1);
    ArrayDimensions=ArrayDimensions(12:end);
    for x=(1:length(Data))
    %This loops put data in a useable form
         if x==2
            AxialDisplacement_C=Data{1,x};
            AxialDisplacement_C=AxialDisplacement_C(12:end);
            AxialDisplacement_C=strrep(AxialDisplacement_C,',','.');
            ArrayDimensions=str2double(AxialDisplacement_C);
            AxialDisplacement=ArrayDimensions;
            ArrayDimensions=zeros(length(Data{1,1}),1);
            ArrayDimensions=ArrayDimensions(12:end);
        end
        if x==3
            AxialForce_C=Data{1,x};
            AxialForce_C=AxialForce_C(12:end);
            AxialForce_C=strrep(AxialForce_C,',','.');
            ArrayDimensions=str2double(AxialForce_C);
            AxialForce=ArrayDimensions;
            ArrayDimensions=zeros(length(Data{1,1}),1);
            ArrayDimensions=ArrayDimensions(12:end);
        end
    end 
    AxialForce=transpose(AxialForce);
    DataCell{i,4}=AxialForce;
end


%% Relate stress to strain using time stamps
for i=(1:length(NameArray))
    
AxialForce=DataCell{i,4};
StrainLine=StrainCell{1,i};
FrameTime=(1/FPS);
SampleTime=(1/MTSHz);

TimeArrayCamera=(ImageNo.*FrameTime)-(FrameTime);
AxialForceTime=(linspace(1,(length(AxialForce)*SampleTime),length(AxialForce))-1);
DifferenceArray=zeros(1,length(AxialForceTime));
AxialForceCut=zeros(1,length(TimeArrayCamera));
NewAxialForceTime=zeros(1,length(TimeArrayCamera));

for j=(1:length(TimeArrayCamera))
    TimeCam=TimeArrayCamera(j);
    for k=(1:length(AxialForceTime))
        Diff=abs(TimeCam-AxialForceTime(k));
        DifferenceArray(k)=Diff;
        if k==1
            Index=k;
        else
            if DifferenceArray(k)<DifferenceArray(k-1)
                Index=k;
            end
        end
    end
    AxialForceCut(j)=AxialForce(Index);
    NewAxialForceTime(j)=AxialForceTime(Index);
end


figure(10)
yyaxis left
plot(TimeArrayCamera,StrainLine,'*')
hold on
yyaxis right
plot(NewAxialForceTime,AxialForceCut)
hold off




%% Get Stress from Force


Area=((Diameters(i)/2)^2)*pi;
Stress=AxialForceCut./Area;


% Get Ultimate Tensile Strength
UltimateTensileStrengthArray(i)=max(Stress);

% Don't use very low stress values to calculate Young's Modulus 
StressCounter=0;
for j=(1:length(Stress)/10)
    temp=Stress(j);
    if temp<StressThreshold
        StressCounter=StressCounter+1;
    end
end
% Store Stress
StressStore{1,i}=Stress;
StressStore{2,i}=StressCounter;
%% Stress vs Strain

figure(FigureCounter)
plot(StrainLine,Stress);
xlabel('Strain')
ylabel('Stress (MPa)')
end

%% Young's Modulus

% Storage arrays
YoungsModulus=zeros(1,length(NameArray));
Resolution=1/InitialLength;

% Find Young's Modulus for each set of data
for i=(1:length(NameArray))
    
FigureCounter=FigureCounter+1;
Strain=StrainCell{1,i};
Stress=StressStore{1,i};
StressCounter=StressStore{2,i};
MarkerStart=0;
MarkerEnd=0;

if AdjustModulusPoints==1
    TempTemp=Strain(1);
StrainStart=StrainStart+TempTemp;
StrainEnd=StrainEnd+TempTemp;
end


% Find initial strain point
for j=(1:length(Strain))
    Temp=Strain(j);
    if Temp>=StrainStart && MarkerStart==0
        MarkerStart=1;
        StressPointS=Stress(j+StressCounter);
        StrainPointS=Strain(j+StressCounter);
        CountStart=j+StressCounter;
    end    
end

% Find end strain point
for j=(1:length(Strain))
    Temp=Strain(j);
    if Temp>=StrainEnd && MarkerEnd==0
        MarkerEnd=1;
        StressPointE=Stress(j+StressCounter);
        StrainPointE=Strain(j+StressCounter);
        CountEnd=j+StressCounter;
    end    
end

if RawData==1
    ArrayRaw1=zeros(1,length(Strain));
    ArrayRaw2=zeros(1,length(Strain));
    CounterRaw1=1;
    CounterRaw2=1;
    for j=(1:length(Strain))
        StrainVal=Strain(j);
        if StrainVal==Resolution
            ArrayRaw1(CounterRaw1)=j;
            if CounterRaw1~=1
                if ArrayRaw1(CounterRaw1)>ArrayRaw1(CounterRaw1-1)
                    IndexMax1=j;
                end 
            end
            CounterRaw1=CounterRaw1+1;
        end
        if StrainVal==Resolution*2
            ArrayRaw2(CounterRaw2)=j;
            if CounterRaw2~=1
                if ArrayRaw2(CounterRaw2)>ArrayRaw2(CounterRaw2-1)
                    IndexMax2=j;
                end 
            end 
            CounterRaw2=CounterRaw2+1;
        end 
    end
    
        StrainPointS=Strain(IndexMax1);
        StrainPointE=Strain(IndexMax2);
        StressPointS=Stress(IndexMax1);
        StressPointE=Stress(IndexMax2);
end

% Find Young's Modulus using gradient of line
Gradient_YM=(StressPointE-StressPointS)/(StrainPointE-StrainPointS);
YoungsModulus(i)=Gradient_YM;

% Plot Graph of Stress vs Strain with modulus line
 limity=(max(Stress))+(0.05*max(Stress));
 limitx=max(Strain);
 MaxTensileStrength=max(Stress);
[Xco1,Yco1]=find(Stress==MaxTensileStrength);

 YcoEnd=Xco1*Gradient_YM;
 
 figure(FigureCounter)
    plot(Strain,Stress)
    title(['Sample ' num2str(i) ])
    xlabel('Strain ')
    ylabel('Stress (MPa)')
    xlim([0,limitx])
    ylim([0,limity])
    hold on
    
    %Plot Young's Modulus
    
    plot([StrainPointS,Xco1],[StressPointS,YcoEnd])
    legend('Stress vs Strain','Young''s Modulus')
    hold off

end
%% Calculations and Plotting of Tensile Strength
%Find average Maximum Tensile Strength and Young's modulus along with their
%standard deviations.
Samples=length(NameArray);
Avg_YM=mean(YoungsModulus);
SD_YM=std(YoungsModulus);

Avg_UTS=mean(UltimateTensileStrengthArray);
SD_UTS=std(UltimateTensileStrengthArray);


% Display values of Maximum Tensile Strength and Young's modulus along with their
%standard deviations

formatprint_YM='Average Young''s Modulus from %3.0f samples is %4.3f +/- %4.3f MPa\n';

fprintf(formatprint_YM,Samples,Avg_YM,SD_YM)

formatprint_UTS='Average Ultimate Tensile Strength from %3.0f samples is %4.3f +/- %4.3f MPa\n';

fprintf(formatprint_UTS,Samples,Avg_UTS,SD_UTS)

%Version15
