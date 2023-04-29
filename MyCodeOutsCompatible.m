%%%                                                     DiazLab Microfluidics Laboratory
%                                                               Malcom Díaz García
%                                                                 30/sept/2021
%                                                    University of Puerto Rico at Mayaguez 
%                                                             HeLa, MDA-MB-231, MCF12A


%                 Main Functions
%1. Import & Save Individual GAMRY uPore Data with or without Cells and Arrange Into Tables
%2. Create .csv Files of Arranged Data to Feed Into R Studio


 %        IMPORTANT: CODE ASSUMES THE FOLLOWING:
 
 %1) User has a folder created for each cell type. Note that
 %     its name must be exact to the one in the code.
 %2)Variables corresponding to addresses must be adapted to the user's PC. Note: Matlab might not be able 
 %     to detect files inside a onedrive folder. In this case move folder to HDD or SSD
 %3)Code is formatted for Gamry equipment, where (1) excel file harbors (1) sheet of data for only (1) cell.
 %     To use with Palmsens4 the loop must be modified.
 %3)To specify data is for no cells the strings must be edited to label it correctly. 
 %4)Function 2 assumes baseline was substracted already. If not variable name must be edited in code manually
 %5)Function 2 assumes baseline substraction was done using code Cell_Outlier_ID_Code. Make sure files are named as per its output name format
                                            
%             Important Variables                                          
%typecounter= counter created to iterate through cell types
%vcounter= counter created to iterate through voltages 
%cscounter= counter created to iterate through cell samples
%dataloc= location of data within PC
%data= data retrieved from Excel

clear all
close all
clc
format 'longg'                                                              %Set Matlab number display format
Vol=[30 60 90 120 150];                                                     %Evaluated Voltages
Cell=["HELA" "MDA" "MCF"];                                                  %List of Cell Types
FLDRPATH = 'D:\Data Review Tesis\';                                         %String for reading Excel document


u=input('Select the function to execute: \n 1. Import & Save Individual GAMRY uPore Data with Cells and Arrange Into Tables \n 2. Create .csv Files of Arranged Data to Feed Into R Studio \n', 's');

if u=="1" 
    y="yes";
%%________________________________________________________________ [1] ______________________________________________________________________________________________________________________________________________________________________________
while y=="yes"
    
SAVESTR='D:\Data Review Tesis\';                                           %PATH TO OUTPUT FILE FOLDER
DATASRCSTR = '\Excel data\P_EIS_';                                         %String for reading Excel document
VUNITSTR= 'MV_ (';                                                         %String for reading Excel document
ZphzR=[]; %Zphz is predefined
ZmodR=[]; %Zmod is predifined 

v=input('Input voltage to extract cell data from under brackets (30mV, 60mV, 90mV, 120mV and 150mV) \n (Hit enter for all) \n');  %Prompt user to enter voltage
              
if isempty(v)                                                              %Enables hitting enter to define 
    v=Vol;                                                                 %all voltages
end
c=input('Input cell types to extract cell data from under brackets (HELA, MDA-MB-231 MCF12A) \n (For MDA-MB-231 and MCF12A type MDA and MCF, respectively) (Hit enter for all) \n', 's'); %Prompt user to enter cell type
if isempty(c)                                                              %Enables hitting enter to define 
    c=Cell;                                                                %all cell types
end

w=input('Was data ran with or without cells (type Wc/Woc)? \n', 's');      %User specifies Wc or Woc and address string is adjusted accordingly
if w=="Wc"
    CELLSTR= '_CELL_';
elseif w=="Woc"
    CELLSTR= '_NO_CELL_'; 
end

n=input('Enter cell sample size per type&voltage \n');                     %Prompt user to enter total number of cells per mV and type uPores evaluated
C=string(split(c)); C=C.';                                                 %Converts c to a string array of numel(c) elements
[Vi,ia]=intersect(Vol, v); vi=ia; vi=vi.';                                 %creates horizontal array vi composed of input voltages' indexes in array Vol
[Ci,ia]=intersect(Cell,C); ci=ia; ci=ci.';                                 %%creates horizontal array ci composed of input cell types' indexes in array Cell

Zmod=[];                                                                   %Creating & Defining Impedance Magnitude as a Matrix
Zphz=[];                                                                   %Creating & DefImpedance Phase Angle as a Matrix
Freq=[];                                                                   %Creating & Defining Frequency as a Matrix

freqloc=strjoin([FLDRPATH, Cell(1), DATASRCSTR, Cell(1), CELLSTR, Vol(1), VUNITSTR, num2str(1), '.xlsx'], '');
Freq=[xlsread(freqloc, 'D71:D126')]; Freq=Freq.';                          %Arranges 1 row frequency matrix
freqfolder=[FLDRPATH 'Frequency'];
mkdir(freqfolder);                                                         %Create Folder to store freq
save([freqfolder '\Freq_.mat'],'Freq');                                    %Store freq in folder
disp('Frequency has been arranged and saved')

for typecounter=1:numel(ci)                                                %cell type loop to change from one to the other
    for vcounter=1:numel(vi)                                               %tested voltage(s) loop to change from one to the other
        for cscounter=1:n                                                  %cell sample loop to stack their data into a table 
            
            dataloc=[FLDRPATH, num2str(Cell(ci(typecounter))), DATASRCSTR, num2str(Cell(ci(typecounter))), CELLSTR, num2str(Vol(vi(vcounter))), VUNITSTR, num2str(cscounter), '.xlsx'];
            datamod=xlsread(dataloc, 'H71:H126');
            dataphz=xlsread(dataloc, 'I71:I126');
            Zmod=[Zmod datamod];
            Zphz=[Zphz dataphz];
        end
                                                                            %SAVE Zmod for all the collected n samples of same mV & type into a .mat table and variable
        qprompt=['Do you want to create a folder to store the ', num2str(Vol(vi(vcounter))), 'mV ', num2str(Cell(ci(typecounter))), ' data in? (yes/no) \n']; %asks user if he want to create a folder to store .mat data in
question=input(qprompt, 's');

if question=="yes"
    savefolder=[SAVESTR num2str(Cell(ci(typecounter))) '\' num2str(Vol(vi(vcounter))) 'mV' '\Matlab Tables'];
    mkdir(savefolder);
    save([savefolder '\Zmod_' w '_' num2str(Cell(ci(typecounter))) '_' num2str(Vol(vi(vcounter))) '.mat'],'Zmod')
    disp('Zmod has been arranged and saved')
    save([savefolder '\Zphz_' w '_' num2str(Cell(ci(typecounter))) '_' num2str(Vol(vi(vcounter))) '.mat'],'Zphz')
    disp('Zphz has been arranged and saved')
else                                                                       %if user types anything other than "yes", .mat tables are stored inside celltype folder 
    save([SAVESTR '\Zmod_' w '_' num2str(Cell(ci(typecounter))) '_' num2str(Vol(vi(vcounter))) '.mat'],'Zmod')  
    disp('Zmod has been arranged and saved')
    save([SAVESTR '\Zphz_' w '_' num2str(Cell(ci(typecounter))) '_' num2str(Vol(vi(vcounter))) '.mat'],'Zphz')
    disp('Zphz has been arranged and saved')
end
        Zmod=[];                                                           %Cleans Zmod and Zphz arrays 
        Zphz=[];                                                           %to accomodate data of other mV/type for the next iteration
    end
end
disp('All the requested data has finished arrangement and saved into PC')

y=input('Prepare another extraction? (yes/no) \n', 's');
end
if y=="no"
    disp('Program ended');
end
%________________________________________________________________ END 1 ______________________________________________________________________________________________________________________________________________________________________ 

elseif u=="2"
%_________________________________________________________________ [2] ______________________________________________________________________________________________________________________________________________________________________________

y="yes";
%loading Frequency file and storing into new variable Freq
load([FLDRPATH 'Frequency\Freq_.mat']); 


while y=="yes"
  
    %User to input all the cell types to compare 
    c=input('Enter cells to compare (HELA MDA MCF) \n', 's');
    C=string(split(c)); C=C.'; %Converts c to a string array of #input elements
    [Ci,ci]=intersect(Cell,C); ci=ci.';


    %ask user to input amount of cells per ctype
    for x=1:numel(ci)
         ask=strjoin(['Enter the amount of ' Cell(ci(x)) ' cells under study \n']);
         cn(x)=input(ask);
    end

    %User selects the component(s) for which the data is to be merged
    comp=input('Select the component to load and compare with frequency \n (Zmod, Zphz or both) \n', 's');
    %input voltage
    v=input('Enter the cell voltage (30mV, 60mV, 90mV, 120mV or 150mV) \n');
    [Vi,vi]=intersect(Vol,v); vi=vi.'; 
    %define variable for the lowest num in the cells to compare // cn will be the num of each cell that will be merged in a mix
    cn=min(cn); 
    %Code for creating cell id column !!!!!!%%%%%
    Cid={'Cells/Freq'};
    for count=1:2*cn
        add=(count);
    Cid=[Cid; add];   
    end
    w=input('Was data ran with cells? (use Wc or Woc) \n', 's');
    askouts=input("Are both cell's data filtered for outliers? \n", 's');
    
    if askouts == "yes"
    askbase=input("Are both cell's data  without baseline? (yes/no) \n", 's'); %asks if data was filtered because if it is, the line to find files must inlcude the specification
    ZphzRb=[]; %Zphz with baseline only is predefined
    ZphzRnob=[]; %Zphz Wo Outliers and Baseline is predefined
    ZmodRb=[]; %Zmod with baseline only is predifined 
    ZmodRnob=[]; %Zmod Wo Outs & Baseline is predefined 
    %Based on chosen component(s), the data will be loaded and merged ---------------------------------------------------------------------
    if comp== "Zphz" %---------------------------------------------------
        
            %Loop for merging only Zphz of selected cells 
            for n=1:numel(ci)
                phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "Nouts25R_Zphz_", w, "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline
                if askbase == "yes"
                    phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZphzNobouts_", num2str(Vol(vi)), '.mat'], ''); %no base
                end
                
                load(phzmatloc);
                if askbase == "yes"                                        %arranges corresponding matrix 
                    ZphzRnob=[ZphzRnob ZphzNObo(:,1:cn)];                  %nob= no baseline nor outliers
                elseif askbase == "no"
                    ZphzRb=[ZphzRb Zphznouts(:,1:cn)];                      %nb= with baseline
                end
            end
        
            %Freq and cell id columns are added 
            if askbase == "yes"
                ZphzRnob=ZphzRnob'; ZphzRnob=num2cell([Freq; ZphzRnob]); ZphzRnob=[Cid ZphzRnob]; %nob= no baseline nor outliers
            elseif askbase == "no"
            ZphzRb=ZphzRb'; ZphzRb=num2cell([Freq; ZphzRb]); ZphzRb=[Cid ZphzRb]; %with base
            end
    
            %Create csv name
            csvname=[];
            for x=1:numel(ci)
            add=[Cell(x) '_'];
            csvname=[csvname add];
            end
        
            %CSV is saved into PC, in a created folder
            adfldr=strjoin([FLDRPATH Cell(ci) 'Mix'],'');                  %name folder 
            mkdir(adfldr);                                                 %create folder
            if askbase == "yes"
               csvloc=strjoin([adfldr '\ZphzRnob_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
               writecell(ZphzRnob, csvloc);
               disp('Zphz Wo Outliers csv has been saved into PC');
            elseif askbase == "no"
                csvloc=strjoin([adfldr '\ZphzRb_EIS_' csvname num2str(Vol(vi)) '.csv'], '');
                writecell(ZphzRb, csvloc);
                disp('Zphz csv has been saved into PC');
            end 
                        
            
             
             
    elseif comp== "Zmod" % ----------------------------------------------
        
            %Loop for merging only Zmod of selected cells 
            for n=1:numel(ci)
                modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "Nouts25R_Zmod_", w, "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline
                if askbase == "yes"
                    modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZmodNobouts_", num2str(Vol(vi)), '.mat'], ''); 
                end
            end 
            load(modmatloc);
            if askbase == "yes"
                ZmodRnob=[ZmodRnob ZmodNObo(:,1:cn)]; %Wo baseline
            elseif askbase == "no"
                ZmodRb=[ZmodRb Zmodnouts(:,1:cn)]; %with baseline
            end
        
            %Freq and cell id columns are added 
             if askbase == "yes"
                ZmodRnob=ZmodRnob'; ZmodRnob=num2cell([Freq; ZmodRnob]); ZmodRnob=[Cid ZmodRnob];
            elseif askbase == "no"
            ZmodRb=ZmodRb'; ZmodRb=num2cell([Freq; ZmodRb]); ZmodRb=[Cid ZmodRb]; 
            end
            
    
            %Create csv name
            csvname=[];
            for x=1:numel(ci)
            add=[Cell(x) '_'];
            csvname=[csvname add];
            end
            
            %CSV is saved into PC. Code creates a folder
            adfldr=strjoin([FLDRPATH Cell(ci) 'Mix'],'');                  %name folder 
            mkdir(adfldr);                                                 %create folder
            if askbase == "yes"
               csvloc=strjoin([adfldr '\ZmodRnob_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
               writecell(ZmodRnob, csvloc);
               disp('Zmod csv created and saved into PC');
            elseif askbase == "no"
                   csvloc=strjoin([adfldr '\ZmodRb_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
                   writecell(ZmodRb, csvloc);
                   disp('Zmod csv created and saved into PC');
            end 
                        
            
            
            
    elseif comp== "both" %-------------------------------------------------
        
        %Loop for merging Zphz and Zmod of selected cells 
            for n=1:numel(ci)
                phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "Nouts25R_Zphz_", w, "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline
                modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "Nouts25R_Zmod_", w, "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline
                if askbase == "yes"
                    phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZphzNobouts_", num2str(Vol(vi)), '.mat'], '');
                    modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZmodNobouts_", num2str(Vol(vi)), '.mat'], '');
                end
            
                load(phzmatloc);
                load(modmatloc);
                if askbase == "yes"
                    ZphzRnob=[ZphzRnob ZphzNObo(:,1:cn)];                            %stack Zmods into a matrix
                    ZmodRnob=[ZmodRnob ZmodNObo(:,1:cn)];                            %stack Zphzs into a matrix
                elseif askbase == "no"
                    ZphzRb=[ZphzRb Zphznouts(:,1:cn)];                            %stack Zmods into a matrix w/ baseline
                    ZmodRb=[ZmodRb Zmodnouts(:,1:cn)];                            %stack Zphzs into a matrix w/ baseline
                end
            end
        
            %Freq and cell id columns are added 
            if askbase == "yes"
                ZphzRnob=ZphzRnob'; ZphzRnob=num2cell([Freq; ZphzRnob]); ZphzRnob=[Cid ZphzRnob];
                ZmodRnob=ZmodRnob'; ZmodRnob=num2cell([Freq; ZmodRnob]); ZmodRnob=[Cid ZmodRnob];
            elseif askbase == "no"
                ZphzRb=ZphzRb'; ZphzRb=num2cell([Freq; ZphzRb]); ZphzRb=[Cid ZphzRb];
                ZmodRb=ZmodRb'; ZmodRb=num2cell([Freq; ZmodRb]); ZmodRb=[Cid ZmodRb];
            end
            
            %Create csv name
            csvname=[];
            for x=1:numel(ci)
            add=[Cell(ci(x)) '_'];
            csvname=[csvname add];
            end
            
            %CSV is saved into PC, in a created folder
            adfldr=strjoin([FLDRPATH, Cell(ci), 'Mix'], '');
            mkdir(adfldr);

            if askbase== "yes"
                modcsvloc=strjoin([adfldr '\ZmodRnob_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
                phzcsvloc=strjoin([adfldr '\ZphzRnob_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
                writecell(ZphzRnob, phzcsvloc);
                disp('Zphz csv created and saved into PC');
                writecell(ZmodRnob, modcsvloc);
                disp('Zmod csv created and saved into PC');
            elseif askbase == "no"
                phzcsvloc=strjoin([adfldr '\ZphzRb_EIS_' csvname num2str(Vol(vi)) '.csv'], '');
                modcsvloc=strjoin([adfldr '\ZmodRb_EIS_' csvname num2str(Vol(vi)) '.csv'], '');
                writecell(ZphzRb, phzcsvloc);
                disp('Zphz csv created and saved into PC');
                writecell(ZmodRb, modcsvloc);
                disp('Zmod csv created and saved into PC');
            end
            
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~if data is NOT filtered for outliers~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif askouts == "no"
    askbase=input("Are both cell's data  without baseline? (yes/no) \n", 's'); %asks if data was filtered because if it is, the line to find files must inlcude the specification
    ZphzRbo=[]; %Zphz with baseline & outs is predefined
    ZphzRnb=[]; %Zphz wo Baseline with outs is predefined
    ZmodRbo=[]; %Zmod with baseline & outs is predifined 
    ZmodRnb=[]; %Zmod wo Baseline with outs is predefined 
    %Based on chosen component(s), the data will be loaded and merged ---------------------------------------------------------------------
    if comp== "Zphz" %---------------------------------------------------
        
            %Loop for merging only Zphz of selected cells 
            for n=1:numel(ci)
                phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", "Zphz_", w, "_", Cell(ci(n)), "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline and outliers
                if askbase == "yes"
                    phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZphzNobase_", num2str(Vol(vi)), '.mat'], ''); %no base with outs
                end
                
                load(phzmatloc);
                if askbase == "yes"                                        %arranges corresponding matrix 
                    ZphzRnb=[ZphzRnb ZphzNOb(:,1:cn)];                     %nb= no baseline
                elseif askbase == "no"
                    ZphzRbo=[ZphzRbo Zphz(:,1:cn)];                          %nb= with baseline
                end
            end
        
            %Freq and cell id columns are added 
            if askbase == "yes"
                ZphzRnb=ZphzRnb'; ZphzRnb=num2cell([Freq; ZphzRnb]); ZphzRnb=[Cid ZphzRnb]; %nb= no baseline 
            elseif askbase == "no"
                ZphzRbo=ZphzRbo'; ZphzRbo=num2cell([Freq; ZphzRbo]); ZphzRbo=[Cid ZphzRbo]; %with base
            end
    
            %Create csv name
            csvname=[];
            for x=1:numel(ci)
            add=[Cell(x) '_'];
            csvname=[csvname add];
            end
        
            %CSV is saved into PC, in a created folder
            adfldr=strjoin([FLDRPATH Cell(ci) 'Mix'],'');                  %name folder 
            mkdir(adfldr);                                                 %create folder
            if askbase == "yes"
               csvloc=strjoin([adfldr '\ZphzRnb_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
               writecell(ZphzRnb, csvloc);
               disp('Zphz W Outliers wo baseline csv has been saved into PC');
            elseif askbase == "no"
                csvloc=strjoin([adfldr '\ZphzRbo_EIS_' csvname num2str(Vol(vi)) '.csv'], '');
                writecell(ZphzRbo, csvloc);
                disp('Zphz W Outliers & Baseline csv has been saved into PC');
            end 
                        
            
             
             
    elseif comp== "Zmod" % ----------------------------------------------
        
            %Loop for merging only Zmod of selected cells 
            for n=1:numel(ci)
                modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", "Zmod_", w, "_", Cell(ci(n)), "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline and outliers
                if askbase == "yes"
                   modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZmodNobase_", num2str(Vol(vi)), '.mat'], ''); %no base with outs 
                end
            end 
            
            load(modmatloc);
            if askbase == "yes"
                ZmodRnob=[ZmodRnb ZmodNOb(:,1:cn)]; %Wo baseline
            elseif askbase == "no"
                ZmodRbo=[ZmodRbo Zmod(:,1:cn)]; %with baseline
            end
        
            %Freq and cell id columns are added 
             if askbase == "yes"
                ZmodRnb=ZmodRnb'; ZmodRnb=num2cell([Freq; ZmodRnb]); ZmodRnb=[Cid ZmodRnb];
            elseif askbase == "no"
            ZmodRbo=ZmodRbo'; ZmodRbo=num2cell([Freq; ZmodRbo]); ZmodRbo=[Cid ZmodRbo]; 
            end
            
    
            %Create csv name
            csvname=[];
            for x=1:numel(ci)
            add=[Cell(x) '_'];
            csvname=[csvname add];
            end
            
            %CSV is saved into PC. Code creates a folder
            adfldr=strjoin([FLDRPATH Cell(ci) 'Mix'],'');                  %name folder 
            mkdir(adfldr);                                                 %create folder
            if askbase == "yes"
               csvloc=strjoin([adfldr '\ZmodRnb_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
               writecell(ZmodRnb, csvloc);
               disp('Zmod wo Baseline w Outliers csv created and saved into PC');
            elseif askbase == "no"
                   csvloc=strjoin([adfldr '\ZmodRbo_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
                   writecell(ZmodRbo, csvloc);
                   disp('Zmod w Baseline w outliers csv created and saved into PC');
            end 
                        
            
            
            
    elseif comp== "both" %-------------------------------------------------
        
        %Loop for merging Zphz and Zmod of selected cells 
            for n=1:numel(ci)
                phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", "Zphz_", w, "_", Cell(ci(n)), "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline and outliers
                modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", "Zmod_", w, "_", Cell(ci(n)), "_", num2str(Vol(vi)), '.mat'], ''); %w/ baseline and outliers
                if askbase == "yes"
                    phzmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZphzNobase_", num2str(Vol(vi)), '.mat'], ''); %no base with outs
                    modmatloc=strjoin([FLDRPATH, Cell(ci(n)), '\', num2str(Vol(vi)), "mV\Matlab Tables\", Cell(ci(n)), "_ZmodNobase_", num2str(Vol(vi)), '.mat'], ''); %no base with outs
                end
            
                load(phzmatloc);
                load(modmatloc);
                if askbase == "yes"
                    ZphzRnb=[ZphzRnb ZphzNOb(:,1:cn)];                            %stack Zmods into a matrix
                    ZmodRnb=[ZmodRnb ZmodNOb(:,1:cn)];                            %stack Zphzs into a matrix
                elseif askbase == "no"
                    ZphzRbo=[ZphzRbo Zphz(:,1:cn)];                            %stack Zmods into a matrix w/ baseline
                    ZmodRbo=[ZmodRbo Zmod(:,1:cn)];                            %stack Zphzs into a matrix w/ baseline
                end
            end
        
            %Freq and cell id columns are added 
            if askbase == "yes"
                ZphzRnb=ZphzRnb'; ZphzRnb=num2cell([Freq; ZphzRnb]); ZphzRnb=[Cid ZphzRnb];
                ZmodRnb=ZmodRnb'; ZmodRnb=num2cell([Freq; ZmodRnb]); ZmodRnb=[Cid ZmodRnb];
            elseif askbase == "no"
                ZphzRbo=ZphzRbo'; ZphzRbo=num2cell([Freq; ZphzRbo]); ZphzRbo=[Cid ZphzRbo];
                ZmodRbo=ZmodRbo'; ZmodRbo=num2cell([Freq; ZmodRbo]); ZmodRbo=[Cid ZmodRbo];
            end
            
            %Create csv name
            csvname=[];
            for x=1:numel(ci)
            add=[Cell(ci(x)) '_'];
            csvname=[csvname add];
            end
            
            %CSV is saved into PC, in a created folder
            adfldr=strjoin([FLDRPATH, Cell(ci), 'Mix'], '');
            mkdir(adfldr);

            if askbase== "yes"
                modcsvloc=strjoin([adfldr '\ZmodRnb_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
                phzcsvloc=strjoin([adfldr '\ZphzRnb_EIS_' csvname num2str(Vol(vi)) '.csv'],'');
                writecell(ZphzRnb, phzcsvloc);
                disp('Zphz csv created and saved into PC');
                writecell(ZmodRnb, modcsvloc);
                disp('Zmod csv created and saved into PC');
            elseif askbase == "no"
                phzcsvloc=strjoin([adfldr '\ZphzRbo_EIS_' csvname num2str(Vol(vi)) '.csv'], '');
                modcsvloc=strjoin([adfldr '\ZmodRbo_EIS_' csvname num2str(Vol(vi)) '.csv'], '');
                writecell(ZphzRbo, phzcsvloc);
                disp('Zphz csv created and saved into PC');
                writecell(ZmodRbo, modcsvloc);
                disp('Zmod csv created and saved into PC');
            end
    end
y=input('Prepare another mix? (yes/no) \n', 's');

end
if y=="no"
    disp('Program ended'); 
end
%_____________________________________________________________ END 2 ______________________________________________________________________________________________________________________________________________________________________ 
end
end

