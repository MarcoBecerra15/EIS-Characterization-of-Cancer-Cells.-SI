%%               Biosensing and Microfluidics Research Laboratory
%                        DiazLab Microfluidics Laboratory
%                                Malcom Díaz García
%                                    12/oct/2021
%                          University of Puerto Rico at Mayaguez
%         Task: Plot Zphz and Zmod vs Frequency uPore Data with and without cells
%          
%               HeLa, MDA-MB-231, MCF12A
%% Set Code
clearvars
close all
clc
format 'longg'

Vol=[30 60 90 120 150];
Cell=["HELA" "MDA" "MCF"];
FLDRPATH = 'D:\Data Review Tesis\';
%% Input desired parameters (cell type, voltage, data with/without cells)

c=input('Enter cells to plot data of (HELA MDA MCF): \n (Hit enter for all) \n', 's');
if isempty(c)                                                               %Enables hitting enter to define 
    c=Cell;                                                                 %all cell types
end
C=string(split(c)); %C=C.'; %Converts c to a string array of #input elements
[Ci,ia]=intersect(Cell,C); ci=ia; ci=ci.';                                 %%creates horizontal array ci composed of input cell types' indexes in array Cell

v=input('Input the tested voltage to plot data for (30mV, 60mV, 90mV, 120mV and 150mV) \n (Hit enter for all) \n');
if isempty(v)                                                               %Enables hitting enter to define all voltages
    v=Vol;
end
[Vi,ia]=intersect(Vol, v); vi=ia; vi=vi.';                                 %creates horizontal array vi composed of input voltages' indexes in array Vol
w=input('Specify whether data was ran with cells or no cells (Use Wc and Woc, correspondingly): \n', 's'); 



%%  Removing the background noise in the data
[a1,a2]=size(ZmodWHeLa);
ZmodAbsHeLa=[];
ZphzAbsHeLa=[];
for i=1:a2
   ZmodAbsHeLa=[ZmodAbsHeLa (ZmodWHeLa(:,i)-ZmodAvgWoHeLa)];                    %Cell Impedance without device background noise -Absolute cell impedance
   ZphzAbsHeLa=[ZphzAbsHeLa (ZphzWHeLa(:,i)-ZphzAvgWoHeLa)];                    %Cell Phaze without device background noise -Absolute cell phaze
    
end

[a3,a4]=size(ZmodWMDA);
ZmodAbsMDA=[];
ZphzAbsMDA=[];
for i=1:a4
   ZmodAbsMDA=[ZmodAbsMDA (ZmodWMDA(:,i)-ZmodAvgWoMDA)];                        %Cell Impedance without device background noise -Absolute cell impedance
   ZphzAbsMDA=[ZphzAbsMDA (ZphzWMDA(:,i)-ZphzAvgWoMDA)];                        %Cell Phaze without device background noise -Absolute cell phaze
    
end

[a5,a6]=size(ZmodWMCF);
ZmodAbsMCF=[];
ZphzAbsMCF=[];
for i=1:a6
   ZmodAbsMCF=[ZmodAbsMCF (ZmodWMCF(:,i)-ZmodAvgWoMCF)];                        %Cell Impedance without device background noise -Absolute cell impedance
   ZphzAbsMCF=[ZphzAbsMCF (ZphzWMCF(:,i)-ZphzAvgWoMCF)];                        %Cell Phaze without device background noise -Absolute cell phaze
    
end

ZmodAbsAvgHeLa= (mean(ZmodAbsHeLa'))';
ZmodAbsAvgMDA= (mean(ZmodAbsMDA'))';
ZmodAbsAvgMCF= (mean(ZmodAbsMCF'))';
ZphzAbsAvgHeLa= (mean(ZphzAbsHeLa'))';
ZphzAbsAvgMDA= (mean(ZphzAbsMDA'))';
ZphzAbsAvgMCF= (mean(ZphzAbsMCF'))';

%% Iteratively Create Plots based on input values

load([FLDRPATH 'Frequency\Freq_.mat']);                                     %loads Frequency matrix

for ctcount=1:numel(ci)                                                     %cell type loop to change from one to the other
    for vcount=1:numel(vi)                                                  %tested voltage(s) loop to change from one to the other
        phzmatloc=strjoin([FLDRPATH, Cell(ci(ctcount)), '\', num2str(Vol(vi(vcount))), 'mV\Matlab Tables\Zphz_',... %Zphz .mat file location
                  w, '_', Cell(ci(ctcount)) '_', num2str(Vol(vi(vcount))), '.mat'], ''); %in terms of loop parameters
        modmatloc=strjoin([FLDRPATH, Cell(ci(ctcount)), '\', num2str(Vol(vi(vcount))), 'mV\Matlab Tables\Zmod_',... %Zphz .mat file location
                  w, '_', Cell(ci(ctcount)) '_', num2str(Vol(vi(vcount))), '.mat'], ''); %in terms of loop parameters
        load(phzmatloc);
        load(modmatloc);
        
        
        %Create Zmod vs Freq plot in x/y log scale
        TZmodSTR=[Cell(ci(ctcount)) w ' Zmod vs Frequency at ' num2str(Vol(vi(vcount))) 'mV' ];
        figure('NumberTitle','off') 
        loglog(Freq,Zmod)
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
 
          DataTipTemplate.DataTipRows(2).Format = '%f';
        %Saves Zmod Plot  
        SAVESTRmod= strjoin([FLDRPATH Cell(ci(ctcount)) '\' num2str(Vol(vi(vcount))) 'mV\' Cell(ci(ctcount)) 'PlotZmodVsF_' w '_' num2str(Vol(vi(vcount))) 'mV'], '');
        saveas(gcf, SAVESTRmod, 'fig');
        
        %Create Zphz vs Freq plot in x log scale  
        TZphzSTR=[Cell(ci(ctcount)) w ' Zphz vs Frequency at ' num2str(Vol(vi(vcount))) 'mV'];
        figure('NumberTitle','off') 
        semilogx(Freq,Zphz)
        title(TZphzSTR,'FontSize',20);
            xlabel('Frequency (Hz)',...        % letra miu \mu
                'FontName','Arial',...         % tipo de letra
                'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                'FontSize',20);                % Tamaño de letra
            ylabel('Angle (Degrees)',...
                'FontName','Arial',...         % tipo de letra
                'FontWeight','b',...           % Normal(n),Light (l),Demi (d),Bold (b) 
                'FontAngle','n',...            % Normal(n),italic (i),oblique(o) 
                'FontSize',20);                % Tamaño de letra
            set(gca, 'fontsize', 20)
            grid on
            
        %Saves Zphz Plot
        SAVESTRphz= strjoin([FLDRPATH Cell(ci(ctcount)) '\' num2str(Vol(vi(vcount))) 'mV\' Cell(ci(ctcount)) 'PlotZphzVsF_' w '_' num2str(Vol(vi(vcount))) 'mV'], '');
        saveas(gcf, SAVESTRphz, 'fig');
        notice=strjoin([Cell(ci(ctcount)) ' plots for' num2str(Vol(vi(vcount))) 'mV' w ' have been saved into PC']);
        disp(notice);
    end
     
end
