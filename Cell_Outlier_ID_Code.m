%%%                                                     DiazLab Microfluidics Laboratory
%                                                               Malcom Díaz García
%                                                                 30/sept/2021
%                                                    University of Puerto Rico at Mayaguez 
%                                                             HeLa, MDA-MB-231, MCF12A
%                                                          Outlier Identification Code

% Main Functions:
% 1. Identify Zmod/Zphz Outliers using Quartiles Method 
% 2. Generate Plots With and Wo Outliers and only Outliers
% 3. Save Data Without both Zmod/Zphz Outliers in .mat format
% 4. Remove baseline from data W/Wo Outliers 

%Note
% 1. Code has "%" in lines for saving c & cpouts, Out Only Plots, Highlight Plots, to save them remove sign and run code
% 2. In 2. enter no when asked if data was filtered 
% 3. Enter no in Function 2. if by 1. there was no outliers in both Zmod/Zphz. Otherwise hit yes as there should be two Nouts files stored 

clear
close
clc
format longg

Vol=[30 60 90 120 150];                                                     %Evaluated Voltages
Cell=["HELA" "MDA" "MCF"];                                                  %List of Cell Types
FLDRPATH = 'D:\Data Review Tesis\';                                         %String for reading Excel document

askfunction=input('Select the function to execute: \n 1. Identify outliers and Generate Data Set without them \n 2. Create Baseline file and Remove from Data Set \n');
%Input desired parameters (cell type, voltage, data with/without cells)

c=input('Enter cells to evaluate for outliers (HELA MDA MCF) (Hit enter for all): \n', 's');
if isempty(c)                                                              %Enables hitting enter to define 
    c=Cell;                                                                %all cell types
end
C=string(split(c)); C=C.';                                                 %Converts c to a string array of numel(c) elements
v=input('Input the tested voltages to plot data for (30mV, 60mV, 90mV, 120mV and 150mV) (Hit enter for all) \n');
if isempty(v)                                                              %Enables hitting enter to define 
    v=Vol;                                                                 %all voltages
end
[Vi,ia]=intersect(Vol, v); vi=ia; vi=vi.';                                 %creates horizontal array vi composed of input voltages' indexes in array Vol
[Ci,ia]=intersect(Cell,C); ci=ia; ci=ci.';                                 %%creates horizontal array ci composed of input cell types' indexes in array Cell


freqloc=strjoin([FLDRPATH, "Frequency\", "Freq_.mat"], ''); 
load(freqloc);

%% 1. Outlier Loop   
if askfunction == 1
   w=input('Specify whether data was ran with cells or no cells (Use Wc and Woc, correspondingly): \n', 's'); 
    for cellcounter=1:numel(ci)                                                %cell type loop to change from one to the other
        for vcounter=1:numel(vi)                                               %tested voltage(s) loop to change from one to the other
    %Zmod outliers

            matloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\Zmod_", w, "_", Cell(ci(cellcounter)), "_", num2str(Vol(vi(vcounter))), '.mat'], ''); %Find Zmod .mat file 
            load(matloc); Zmsample=Zmod(1:21,:); %Zmod with only the descrete frequencies to be evaluated
            n=size(Zmsample,1);                                                % # of frequencies (rows) evaluated
            m=0.25*n;                                                          %Fail Criteria 
            Zotest=isoutlier(Zmsample, 'quartiles', 2); 
            Zscore=sum(Zotest,1);                                              %Sum the # of f's that each cell classified as outlier
            couts=find(Zscore >= m)                                            %Find the cells that were outliers on 50% or more of the total frequencies evaluated
            Zmouts=Zmod(:,couts); %Zmod outliers only array 
            Zmodno=Zmod; Zmodno(:,couts)=[]; %Zmod without zmod outliers array

            %saveCoutsloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter)), "mV\Matlab Tables\", c, "outs", num2str(100*m/n), "R_Zmod_", w, "_", num2str(Vol(vi(vcounter))), '.mat'],''); %save Cout array into matfile
            %save(saveCoutsloc,'couts');


            %Create Zmod Outliers Only vs Freq plot in x/y log scale
            TZmodSTR=[Cell(ci(cellcounter)) w " Zmod Outliers by Quartile Method vs Frequency at " num2str(100*m/n) "% Pass Requirement: " num2str(Vol(vi(vcounter))) 'mV' ];
            figure('NumberTitle','off') 
            loglog(Freq,Zmouts)
            title(TZmodSTR,'FontSize',20);
                    xlabel('Frequency (Hz)',...        % letra miu \mu
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    ylabel('Impedance (Ohms)',...
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    set(gca, 'fontsize', 20)
                    grid on

                    %Saves Zmod Outliers Only Plot if user wants to
             %askSTR=strjoin(["Save Outliers Only Plot of ", num2str(Vol(vi(vcounter))),"mV ", Cell(ci(cellcounter)), " Zmod? (yes/no) \n"], '');
             %ask=input(askSTR,'s');
             %if ask == "yes"
                    %SAVESTRmod= strjoin([FLDRPATH Cell(ci(cellcounter)) "\" num2str(Vol(vi(vcounter))) 'mV\' Cell(ci(cellcounter)) 'PlotR' num2str(100*m/n) 'ZmodoutsVsF_' w '_' num2str(Vol(vi(vcounter))) 'mV'], '');
                    %saveas(gcf, SAVESTRmod, 'fig');
             %end

             %Create Zmod without Outliers vs Freq plot in x/y log scale
             TZmodSTR=[Cell(ci(cellcounter)) w " Zmod Without Outliers vs Frequency at " num2str(100*m/n) "% Pass Requirement: " num2str(Vol(vi(vcounter))) 'mV' ];
             figure('NumberTitle','off') 
             loglog(Freq,Zmodno)
             title(TZmodSTR,'FontSize',20);
                    xlabel('Frequency (Hz)',...        % letra miu \mu
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    ylabel('Impedance (Ohms)',...
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    set(gca, 'fontsize', 20)
                    grid on

                    %Saves Zmod Wo Outliers Plot if user wants to 
             askSTR=strjoin(["Save Plot Without Outliers of ", num2str(Vol(vi(vcounter))),"mV ", Cell(ci(cellcounter)), " Zmod? (yes/no) \n"], '');
             ask=input(askSTR,'s');
             if ask == "yes"
                SAVESTRmod= strjoin([FLDRPATH Cell(ci(cellcounter)) "\" num2str(Vol(vi(vcounter))) 'mV\' Cell(ci(cellcounter)) 'PlotR' num2str(100*m/n) 'ZmodnoVsF_' w '_' num2str(Vol(vi(vcounter))) 'mV'], '');
                saveas(gcf, SAVESTRmod, 'fig');
             end

            %Create plot highlighting the outliers among the data
            TZmodSTR=[Cell(ci(cellcounter)) w " Zmod with Highlighted Outliers (Red) vs Frequency at " num2str(100*m/n) "% Pass Requirement: " num2str(Vol(vi(vcounter))) 'mV' ];
            figure('NumberTitle','off')
            loglog(Freq, Zmodno, 'b', Freq, Zmouts, 'r')
             title(TZmodSTR,'FontSize',20);
                    xlabel('Frequency (Hz)',...        % letra miu \mu
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    ylabel('Impedance (Ohms)',...
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    set(gca, 'fontsize', 20)
                    grid on

                    %Saves Zmod with Highlighted Outliers Plot    
             %askSTR=strjoin(["Save Plot with Highlighted Outliers of ", num2str(Vol(vi(vcounter))),"mV ", Cell(ci(cellcounter)), " Zmod? (yes/no) \n"], '');
             %ask=input(askSTR,'s');
             %if ask == "yes"
                    %SAVESTRmod= strjoin([FLDRPATH Cell(ci(cellcounter)) "\" num2str(Vol(vi(vcounter))) 'mV\' Cell(ci(cellcounter)) 'PlotR' num2str(100*m/n) 'ZmodRedOutsVsF_' w '_' num2str(Vol(vi(vcounter))) 'mV'], '');
                    %saveas(gcf, SAVESTRmod, 'fig'); 
             %end

     %Zphz Outliers---------------------------------------------------------

            matloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\Zphz_", w, "_", Cell(ci(cellcounter)), "_", num2str(Vol(vi(vcounter))), '.mat'], ''); %Find Zphz
            load(matloc); Zpsample=Zphz(1:42,:);                               %Zphz with only the descrete frequencies to be evaluated
            n=size(Zpsample,1);                                                % # of frequencies (rows) evaluated
            m=0.25*n;                                                          %Fail Criteria 
            Zotest=isoutlier(Zpsample, 'quartiles', 2); 
            Zscore=sum(Zotest,1); %Sum the # of f's that each cell classified as outlier
            cpouts=find(Zscore >= m) %Find the cells that were outliers on 50% or more of the total frequencies evaluated
            Zpouts=Zphz(:,cpouts); %Zphz outliers only 
            Zphzno=Zphz; Zphzno(:,cpouts)=[]; %Zphz without outliers

            %saveCoutsloc=strjoin([FLDRPATH, c, '\', v, "mV\Matlab Tables\", c, "outs", num2str(100*m/n), "R_Zphz_", w, "_", v, '.mat'],''); %save Cout array into matfile
            %save(saveCoutsloc,'cpouts');


            %Create Zphz Outliers Only vs Freq plot in x/y log scale
            TZphzSTR=[Cell(ci(cellcounter)) w " Zphz Outliers by Quartile Method vs Frequency at " num2str(100*m/n) "% Pass Requirement: " num2str(Vol(vi(vcounter))) 'mV' ];
            figure('NumberTitle','off') 
            semilogx(Freq,Zpouts)
            title(TZphzSTR,'FontSize',20);
                    xlabel('Frequency (Hz)',...        % letra miu \mu
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    ylabel('Impedance (Ohms)',...
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    set(gca, 'fontsize', 20)
                    grid on


                  %Saves Zphz Outliers Only Plot if user wants to  
            %askSTR=strjoin(["Save Outliers Only Plot of ", num2str(Vol(vi(vcounter))),"mV ", Cell(ci(cellcounter)), " Zphz? (yes/no) \n"], '');
             %ask=input(askSTR,'s');
             %if ask == "yes"
                    %SAVESTRphz= strjoin([FLDRPATH Cell(ci(cellcounter)) "\" num2str(Vol(vi(vcounter))) 'mV\' Cell(ci(cellcounter)) 'PlotR' num2str(100*m/n) 'ZphzoutsVsF_' w '_' num2str(Vol(vi(vcounter))) 'mV'], '');
                    %saveas(gcf, SAVESTRphz, 'fig');
             %end

            %Create Zphz Wo Outliers vs Freq plot in x log scale
            TZphzSTR=[Cell(ci(cellcounter)) w " Zphz Without Outliers vs Frequency at " num2str(100*m/n) "% Pass Requirement: " num2str(Vol(vi(vcounter))) 'mV' ];
            figure('NumberTitle','off') 
            semilogx(Freq,Zphzno)
            title(TZphzSTR,'FontSize',20);
                    xlabel('Frequency (Hz)',...        % letra miu \mu
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    ylabel('Impedance (Ohms)',...
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    set(gca, 'fontsize', 20)
                    grid on


                   %Saves Zmod Wo Outliers Plot if user wants to 
             askSTR=strjoin(["Save Plot Without Outliers of ", num2str(Vol(vi(vcounter))),"mV ", Cell(ci(cellcounter)), " Zphz? (yes/no) \n"], '');
             ask=input(askSTR,'s');
             if ask == "yes"
                SAVESTRphz= strjoin([FLDRPATH Cell(ci(cellcounter)) "\" num2str(Vol(vi(vcounter))) 'mV\' Cell(ci(cellcounter)) 'PlotR' num2str(100*m/n) 'ZphznoVsF_' w '_' num2str(Vol(vi(vcounter))) 'mV'], '');
                saveas(gcf, SAVESTRphz, 'fig');
             end


             %Create plot highlighting the outliers among the data
            TZphzSTR=[Cell(ci(cellcounter)) w " Zphz with Highlighted Outliers (Red) vs Frequency at " num2str(100*m/n) "% Pass Requirement: " num2str(Vol(vi(vcounter))) 'mV' ];
            figure('NumberTitle','off')
            semilogx(Freq, Zphzno, 'b', Freq, Zpouts, 'r')
             title(TZphzSTR,'FontSize',20);
                    xlabel('Frequency (Hz)',...        % letra miu \mu
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    ylabel('Impedance (Ohms)',...
                     'FontName','Arial',...         % tipo de letra
                     'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                     'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                     'FontSize',20);                % Tamaño de letra
                    set(gca, 'fontsize', 20)
                    grid on


            %Saves Zphz with Highlighted Outliers Plot    
             %askSTR=strjoin(["Save Plot with Highlighted Outliers of ", Vol(vi(vcounter)),"mV ", Cell(ci(cellcounter)), " Zphz? (yes/no) \n"], '');
             %ask=input(askSTR,'s');
             %if ask == "yes"
                    %SAVESTRphz= strjoin([FLDRPATH Cell(ci(cellcounter)) "\" Vol(vi(vcounter)) 'mV\' Cell(ci(cellcounter)) 'PlotR' num2str(100*m/n) 'ZphzRedOutsVsF_' w '_' Vol(vi(vcounter)) 'mV'], '');
                    %saveas(gcf, SAVESTRphz, 'fig'); 
             %end


     %Combine Zmod/Zphz Outliers 
            comouts=intersect(couts,cpouts);
            uncomouts=setxor(couts, cpouts);
            zmpouts= sort([comouts uncomouts])
            Zphznouts=Zphz; Zphznouts(:,zmpouts)=[]; %Zphz without outliers from both Zhphz and Zmod (worst case)
            Zmodnouts=Zmod; Zmodnouts(:, zmpouts)=[]; %Zmod without outliers from both Zphz and Zmod (worst case)

            saveNopouts=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "Nouts", num2str(100*m/n), "R_Zphz_", w, "_", num2str(Vol(vi(vcounter))), '.mat'],''); %save filtered Zphz array into matfile
            save(saveNopouts,'Zphznouts');

            saveNomouts=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "Nouts", num2str(100*m/n), "R_Zmod_", w, "_", num2str(Vol(vi(vcounter))), '.mat'],''); %save filtered Zmod array into matfile
            save(saveNomouts,'Zmodnouts');
        end
    end

%% 2. Baseline Generation and Removal 
elseif askfunction == 2

    askouts=input('Is the data filtered for outliers? (yes/no) \n', 's');
    
        for cellcounter=1:numel(ci)                                                %cell type loop to change from one to the other
            for vcounter=1:numel(vi)                                               %tested voltage(s) loop to change from one to the other
                %Calling data
                ncpmatloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\Zphz_Woc_", Cell(ci(cellcounter)), "_", num2str(Vol(vi(vcounter))), '.mat'], '');
                ncmmatloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\Zmod_Woc_", Cell(ci(cellcounter)), "_", num2str(Vol(vi(vcounter))), '.mat'], '');
                phzmatloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\Zphz_Wc_", Cell(ci(cellcounter)), "_", num2str(Vol(vi(vcounter))), '.mat'], ''); %Find Zphz
                modmatloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\Zmod_Wc_", Cell(ci(cellcounter)), "_", num2str(Vol(vi(vcounter))), '.mat'], ''); %Find Zmod .mat file
                if askouts == "yes"
                    phzmatloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "Nouts25R_Zphz_Wc_", num2str(Vol(vi(vcounter))), '.mat'],''); 
                    modmatloc=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "Nouts25R_Zmod_Wc_", num2str(Vol(vi(vcounter))), '.mat'],'');
                end
                
                load(ncpmatloc); load(ncmmatloc);
                ncZmod=Zmod; ncZphz=Zphz;                                  %Since nc and cell containing data may be named by the same variable (Zmod/Zphz), the value of nc is passed unto new variable before being overwritten
                load(phzmatloc); load(modmatloc);                          %Load data containing cells
                askSTR=strjoin(['Enter ' Cell(ci(cellcounter)) ' cells at ' Vol(vi(vcounter)) 'mV to be used to calculate baseline (brackets [] required) \n'], '');
                %Call and make Baseline
                basecells=input(askSTR);                                   %Enter cells to be evaluated for baseline
                pbase=ncZphz(:,basecells);                                  %Store Zphz of basecells                                
                mbase=ncZmod(:,basecells);                                  %Store Zmod of basecells

                pbaseline=mean(pbase,2);                                   %Phz Baseline is the avg of the eval. cells
                mbaseline=mean(mbase,2);                                   %Mod Baseline is the avg of the eval. cells
                y=input('Save Baselines into PC? \n', 's');                   %Optionally save baselines 
                if y== "yes"
                    savembaseSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_Zmod_Baseline_", num2str(Vol(vi(vcounter))), '.mat'],'');
                    savepbaseSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_Zphz_Baseline_", num2str(Vol(vi(vcounter))), '.mat'],'');
                    savembasedataSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_Zmod_BaseCdata_", num2str(Vol(vi(vcounter))), '.mat'],'');
                    savepbasedataSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_Zphz_BaseCdata_", num2str(Vol(vi(vcounter))), '.mat'],'');
                    save(savembaseSTR, 'mbaseline'); save(savepbaseSTR, 'pbaseline'); save(savembasedataSTR, 'mbase'); save(savepbasedataSTR, 'pbase');
                    disp('Baselines have been saved into PC')
                end
                %Substract baseline from the components and save them
                ZmodNOb= Zmod - mbaseline;
                ZphzNOb= Zphz - pbaseline;
                savedatamSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_ZmodNobase_", num2str(Vol(vi(vcounter))), '.mat'],'');
                savedatapSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_ZphzNobase_", num2str(Vol(vi(vcounter))), '.mat'],'');
                if askouts == "yes"
                    ZmodNObo= Zmodnouts - mbaseline;                       %Zmod without baseline nor outliers
                    ZphzNObo= Zphznouts - pbaseline;                      %Zphz without baseline nor outliers
                    savedatamSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_ZmodNobouts_", num2str(Vol(vi(vcounter))), '.mat'],'');
                    savedatapSTR=strjoin([FLDRPATH, Cell(ci(cellcounter)), '\', num2str(Vol(vi(vcounter))), "mV\Matlab Tables\", Cell(ci(cellcounter)), "_ZphzNobouts_", num2str(Vol(vi(vcounter))), '.mat'],'');
                    save(savedatapSTR, 'ZphzNObo'); save(savedatamSTR, 'ZmodNObo'); %Save data filtered into PC
                    disp('Data without Baseline saved into PC');
                else
                    save(savedatapSTR, 'ZphzNOb'); save(savedatamSTR, 'ZmodNOb'); %Save data not filtered into PC
                    disp('Data without Baseline saved into PC');
                end
            end
        end
 
       
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    