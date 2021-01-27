%% mRNA and Protein Accumulation from Spot fluorescence (from the MS2 spot data)
%Adapted from code by YJK
%Note: This will plot spot fluorescence, mRNA accumulation, and predicted
%protein fluorescence for all cells tracked thus far.
%On Nov 12-2019 7 cells have been tracked in two embryos.

%Spot data is saved as follows:
%Each cell has it's own file called Embryo#Cell#.mat
%Within the file is:
    %Fluo - spot intensities for frames for which we have data.
    %ElapsedTime - storing the elapsed time of the recorded frames
    %FluoError - FluoError stored in vector positions corresponding to the Fluo
        %vector. Since two files were sometimes combined there were different 
        %errors for each file.
    


%Set parameters
gammaM = 4.57; %mRNA half-life in min min
lambdam=log(2)./gammaM; %degradation rate in # per min
PProduct= 4.2; %Translation rate Mature mRNA per minute
lambdap=20*lambdam; %Prot deg rate in # per minute
Nt=300; %Total # of time points - 300 is more than enough for the data right now.
dt=0.5; % time steps in minutes
ScaleFac=4; %Scaling factor for plotting later
Protoff=5; %offset to plot protein trace to account for maturation of mRNA,...
            %protein, and fluorophore. In minutes.

%Load appropriate data 
%insert folder locations here if you want.
CellNames=dir('Embryo*Cell*.mat');

%Loop over the cells
for k=1:length(CellNames);
    %Load the ith data
    load(CellNames(k).name);
    
    %Initialize vectors 
    AccumulatedmRNA=zeros(1,Nt);
    AccumulatedmRNAMin=zeros(1,Nt);
    AccumulatedmRNAMax=zeros(1,Nt);
    AccumulatedProtein=zeros(1,Nt);
    AccumulatedProteinMin=zeros(1,Nt);
    AccumulatedProteinMax=zeros(1,Nt);
    PlotTime=zeros(1,Nt);
    PlotmRNA=zeros(1,Nt);
    PlotmRNAMin=zeros(1,Nt);
    PlotmRNAMax=zeros(1,Nt);
    PlotError=zeros(1,Nt);
   
    %Loop over time
    for i=2:Nt;
        
        %add time points
        PlotTime(i)=PlotTime(i-1)+dt;

            %Check if this time point has fluorescent spot data
            if ~isempty(find(ElapsedTime==PlotTime(i)));
                %If there's fluoresence, include it in the spot
                %fluorescence to plot
                PlotmRNA(i)=Fluo(find(ElapsedTime==PlotTime(i)));
                %Get the error for this timepoint
                PlotError(i)=FluoError(find(ElapsedTime==PlotTime(i)));
                %Get error bars 
                PlotmRNAMax(i)=PlotmRNA(i)+PlotError(i);
                PlotmRNAMin(i)=PlotmRNA(i)+PlotError(i); 
                %Make sure that it doesn't go negative - this will mess up
                %the accumulated mRNA and Protein
                    if PlotmRNAMin(i)<0;
                        PlotmRNAMin(i)=0;
                    end
            end
        %Based on number of mRNA in previous time point calc mRNA degraded,
        %mRNA made, and protein made
        AccumulatedmRNA(i)=AccumulatedmRNA(i-1)+PlotmRNA(i)-...
            AccumulatedmRNA(i-1)*lambdam*dt;
        %Calculate min and max values
        AccumulatedmRNAMin(i)=AccumulatedmRNAMin(i-1)+PlotmRNAMin(i)-...
            AccumulatedmRNAMin(i-1)*lambdam*dt;
         AccumulatedmRNAMax(i)=AccumulatedmRNAMax(i-1)+PlotmRNAMax(i)-...
            AccumulatedmRNAMax(i-1)*lambdam*dt;
        %Get erro bar lengths
        AccumulatedmRNAMinError(i)=abs(AccumulatedmRNAMin(i)-AccumulatedmRNA(i));
        AccumulatedmRNAMaxError(i)=AccumulatedmRNAMax(i)-AccumulatedmRNA(i);
        
        
        
        AccumulatedProtein(i)=AccumulatedProtein(i-1)+PProduct*AccumulatedmRNA(i-1)*dt-...
            AccumulatedProtein(i-1)*lambdap*dt;  
        %Calculate min and max values
        AccumulatedProteinMin(i)=AccumulatedProteinMin(i-1)+PProduct*AccumulatedmRNAMin(i-1)*dt-...
            AccumulatedProteinMin(i-1)*lambdap*dt; 
        AccumulatedProteinMax(i)=AccumulatedProteinMax(i-1)+PProduct*AccumulatedmRNAMax(i-1)*dt-...
            AccumulatedProteinMax(i-1)*lambdap*dt; 
        %Get erro bar lengths
        AccumulatedProteinMinError(i)=abs(AccumulatedProteinMin(i)-AccumulatedProtein(i));
        AccumulatedProteinMaxError(i)=AccumulatedProteinMax(i)-AccumulatedProtein(i);
        

    end

    
    %Now plot everything
    figure(k);
    %Transcription spot
    yyaxis left
    ylabel('Fluorescence Intensity (a.u.)');
    errorbar(PlotTime,PlotmRNA*ScaleFac,PlotError*ScaleFac','-r','LineWidth',1);
    hold on
    %Protein
    yyaxis left
    ylabel('Fluorescence Intensity (a.u.)');
    errorbar(PlotTime+Protoff,AccumulatedProtein, AccumulatedProteinMinError,...
         AccumulatedProteinMaxError,'b','LineWidth',1);
    ylim([0 inf])
    hold on
    %Accumulated mRNA
    yyaxis right
    ylabel('mRNA count');
    errorbar(PlotTime,AccumulatedmRNA,AccumulatedmRNAMinError,AccumulatedmRNAMaxError,'-g','LineWidth',1);
    ylim([0 inf])
    hold off

    xlim([0 150])
    xlabel('Time (min)')
    legend('Nascent mRNA','Predicted Protein Signal','Predicted Accumulated mRNA');
    title(CellNames(k).name(1:end-4));
%     StandardFigure(gcf,gca);


end



%% Plot only MS2


%close all
ScaleFac=1/1000;


%7, 2

k=7;
load(CellNames(k).name);

PlotTime=zeros(1,Nt);
PlotmRNA=zeros(1,Nt);
PlotmRNAMin=zeros(1,Nt);
PlotmRNAMax=zeros(1,Nt);
PlotError=zeros(1,Nt);
%Loop over time
for i=2:Nt;

    %add time points
    PlotTime(i)=PlotTime(i-1)+dt;

        %Check if this time point has fluorescent spot data
        if ~isempty(find(ElapsedTime==PlotTime(i)));
            %If there's fluoresence, include it in the spot
            %fluorescence to plot
            PlotmRNA(i)=Fluo(find(ElapsedTime==PlotTime(i)));
            %Get the error for this timepoint
            PlotError(i)=FluoError(find(ElapsedTime==PlotTime(i)));
        end
end


figure(1);
load(CellNames(k).name);
%Transcription spot
ylabel('number of RNAP molecules (a.u.)');
%PlotHandle=errorbar(PlotTime,PlotmRNA*ScaleFac,PlotError*ScaleFac','.-k','LineWidth',1);
%PlotHandle.CapSize=0;
PlotHandle=plot(PlotTime,PlotmRNA*ScaleFac,'-k');
xlim([0 125])
ylim([0,5.5])
xlabel('time (min)')
%title(CellNames(k).name(1:end-4));
set(gca,'YTick',[0:1:6])
StandardFigurePBoC(PlotHandle,gca)






k=2;
load(CellNames(k).name);

PlotTime=zeros(1,Nt);
PlotmRNA=zeros(1,Nt);
PlotmRNAMin=zeros(1,Nt);
PlotmRNAMax=zeros(1,Nt);
PlotError=zeros(1,Nt);
%Loop over time
for i=2:Nt;

    %add time points
    PlotTime(i)=PlotTime(i-1)+dt;

        %Check if this time point has fluorescent spot data
        if ~isempty(find(ElapsedTime==PlotTime(i)));
            %If there's fluoresence, include it in the spot
            %fluorescence to plot
            PlotmRNA(i)=Fluo(find(ElapsedTime==PlotTime(i)));
            %Get the error for this timepoint
            PlotError(i)=FluoError(find(ElapsedTime==PlotTime(i)));
        end
end


figure(2);
load(CellNames(k).name);
%Transcription spot
ylabel('number of RNAP molecules (a.u.)');
%PlotHandle=errorbar(PlotTime,PlotmRNA*ScaleFac,PlotError*ScaleFac','.-k','LineWidth',1);
%PlotHandle.CapSize=0;
PlotHandle=plot(PlotTime,PlotmRNA*ScaleFac,'-k');
xlim([0 125])
ylim([0,5.5])
xlabel('time (min)')
%title(CellNames(k).name(1:end-4));
set(gca,'YTick',[0:1:6])
StandardFigurePBoC(PlotHandle,gca)






%% Browse through MS2 traces

%close all
ScaleFac=1/1000;




for k=1:length(CellNames)

    load(CellNames(k).name);

    PlotTime=zeros(1,Nt);
    PlotmRNA=zeros(1,Nt);
    PlotmRNAMin=zeros(1,Nt);
    PlotmRNAMax=zeros(1,Nt);
    PlotError=zeros(1,Nt);
    %Loop over time
    for i=2:Nt;

        %add time points
        PlotTime(i)=PlotTime(i-1)+dt;

            %Check if this time point has fluorescent spot data
            if ~isempty(find(ElapsedTime==PlotTime(i)));
                %If there's fluoresence, include it in the spot
                %fluorescence to plot
                PlotmRNA(i)=Fluo(find(ElapsedTime==PlotTime(i)));
                %Get the error for this timepoint
                PlotError(i)=FluoError(find(ElapsedTime==PlotTime(i)));
            end
    end
    
    
    figure(k)
    clf
    load(CellNames(k).name);
    %Transcription spot
    ylabel('number of RNAP molecules (a.u.)');
    %PlotHandle=errorbar(PlotTime,PlotmRNA*ScaleFac,PlotError*ScaleFac','.-k','LineWidth',1);
    %PlotHandle.CapSize=0;
    PlotHandle=plot(PlotTime,PlotmRNA*ScaleFac,'.-k');
    xlim([0 125])
    ylim([0,5.5])
    xlabel('time (min)')
    %title(CellNames(k).name(1:end-4));
    set(gca,'YTick',[0:1:6])
    StandardFigure(PlotHandle,gca)
    title(['Trace ', num2str(k)])


end

