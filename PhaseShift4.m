% Make a table with the following columns: residue# TOCSYHB2int TOCSYHB3int
% TOCSYerror HN(CO)HBHB2int HN(CO)HBHB3int HN(CO)HBerror HNHBHB2int
% HNHBHB3int HNHBerror Residue(Single Letter Code) Degenerate(y/n)
% PhaseShift(y/n)
% PhiAngle PsiAngle Characteristic(a/b/av)
% set f1=table prior to running this script
% Minimization variable Int(1), Int(2), Int(3). These are scaling parameters to make predicted peak intensities correlated with experimental.
% Turns out that Int equals the maximum intensity of a peak in each spectrum if one has complete magnetization transfer (sin 90 degrees)
% For Karplus equation parameters see Perez and Schmidt JACS (2001) 123:7081-7093

%Because input table is a double array, non-numbers are not loaded in so
%they have to be replaced with # found below:

%Residue: 
%G -> 1, A -> 2, V -> 3, L -> 4, I -> 5, S -> 6, T -> 7, F -> 8, Y -> 9, H -> 10,
%W -> 11, C -> 12, M -> 13, K -> 14, R -> 15, N -> 16, Q -> 17, D -> 18, E -> 19, P -> 20

%True or false
% 1 for true, 0 for false

%structure prediction: 
%Alpha-helix: 1, beta-strand: 2, Average - 3

%Makes a new directory to store plots
currentDir = pwd;
mkdir RigidPlot
mkdir BetaBranchPlot
mkdir FlexibleDegenPlots
formatDir = '%s/%s';
rigidDir = sprintf(formatDir,currentDir,'RigidPlot');
betaBranchDir = sprintf(formatDir,currentDir,'BetaBranchPlot');
flexibleDegenDir = sprintf(formatDir,currentDir,'FlexibleDegenPlots');

% make empty array to fill with error function minima for contour plot 
[t,z] = meshgrid(1:1:101, 1:1:101); 
L = zeros(101,101); 

TOCSYmix = 0.017;
HNCOHBmix = 0.025;
HNHBmix = 0.038;

%enables script to use phi, psi angles to determine mixed, alpha helix or
%beta sheet characteristic
autoCharacteristic = true;

%Allows for setting plot font characteristics
fontType = 'Arial';
fontSize = 14;

%Will spit out an error 
phiMin = -180;
phiMax = -20;
%Will determine region to average values. Anything below is alpha, anything
%above is beta
psiMin = -10;
psiMax = 50;

counter = 0;
f1=f1(:,:);
nn=length(f1(:,1));
prompt_model = 0;
C0HA2 = 0;

TBetaBranch = zeros(nn,8);
TFlexOrDegen = zeros(nn, 11);
TRigid = zeros(nn,13);
a3 = zeros(360,13);
d = zeros(100,11);
g = zeros(100,8);

prompt_coefficient = input('Enter 1 for Perez coefficients or 2 for Case coefficients: '); 

%Should be nn - 1 to ensure index not out of bounds
while counter <= nn - 1;
    y=f1(counter + 1,:);
    %y(x) - x specifies column of data
    %y(1) - Residue #
    %y(2) - TOCSYHB2int
    %y(3) - TOCSYHB3int
    %y(4) - TOCSYerror
    %y(5) - HN(CO)HBHB2int
    %y(6) - HN(CO)HBHB3int
    %y(7) - HN(CO)HBError
    %y(8) - HNHBHB2int
    %y(9) - HNHBHB3int
    %y(10) - HNHBError
    %y(11) - Residue type
    %y(12) - Degenerate (y/n)
    %y(13) - Phase shift (y/n)
    %y(14) - Phi angle
    %y(15) - Psi angle
    %y(16) -Manual characteristic override
    
    %Counter will be the row of the current, tested residue
    %this is the testing area
    if y(11) == 1 || y(11) == 2 || y(11) == 20;
        %Insert save residue name and number here
        %fprintf('this is an alanine \n')
        %This skips a calculation that does not exist
    else
        if y(12) == 1
            prompt_model = 4;
        elseif y(11) == 3 || y(11) == 5 || y(11) == 7
            prompt_model = 3;
        else
            prompt_model = 1;
        end
        
        %Sets up Perez Coeff. - only have to set this up once
        
        if prompt_coefficient == 1 && C0HA2 == 0           
            C0HA2 = 5.83;
            C0HA3 = C0HA2;
			C1HA2 = -1.37;
            C1HA3 = C1HA2;
			C2HA2 = 3.61;
            C2HA3 = C2HA2;
			dHA2 = 0;
            dHA3 = 0;
            
			C1N2 = -0.75;
            C1N3 = C1N2;
			C2N2 = 1.15;
            C2N3 = C2N2;
			dN2 = 0;
            dN3 = 0;
					  
			C0C2 = 3.32;
            C0C3 = C0C2;
			C1C2 = -1.58;
            C1C3 = C1C2;
			C2C2 = 2.01;
            C2C3 = C2C2;
			dC2 = 0;
            dC3 = 0;
            
            %For beta branched
            
            C0HA = 5.83;
			C1HA = -1.37;
			C2HA = 3.61;
			dHA = 0;
			  
			C0N = 2.22;
			C1N = -0.75;
            C2N = 1.15;
			dN = 0;
			  
			C0C = 3.32;
			C1C = -1.58;
			C2C = 2.01;
			dC = 0;
        elseif prompt_coefficient == 2
            %Sets up Case 
            %Sets up characteristic if enabled
            if autoCharacteristic
                
                if y(14) <= phiMin || y(14) >= phiMax
                    fprintf('Error: phi angle is out of acceptable bounds! \n')
                else
                    if y(15) <= psiMin
                        prompt_avg = 1; %Alpha
                    elseif y(15) >= psiMax
                        prompt_avg = 2; %Beta
                    else
                        prompt_avg = 3; %Av
                    end
                end
            else
                %Sets up average to be input
                prompt_avg = y(16);
            end 
            %Sets up coeff.
            
            if y(11) == 6
                if prompt_avg == 1 %Alpha
                    C0N2 = -2.10;
					C1N2 = 0.45; 
					C2N2 = -2.19; 
					dN2 = -16.50; 
						 		   
					C0N3 = -2.30; 
					C1N3 = 0.67; 
					C2N3 = -2.28; 
					dN3 = 7.60; 
                    
                elseif prompt_avg == 2 %Beta
                    C0N2 = -2.62;
					C1N2 = -0.25; 
					C2N2 = -2.77; 
					dN2 = -14.70; 
						 		   
					C0N3 = -2.77; 
					C1N3 = -0.06; 
					C2N3 = -2.84; 
					dN3 = 11.20;
                    
                else %Mixed
                    C0N2 = -2.36;
                    C1N2 = 0.10; 
                    C2N2 = -2.48; 
                    dN2 = -15.60;

                    C0N3 = -2.54; 
                    C1N3 = 0.31; 
                    C2N3 = -2.56; 
                    dN3 = 9.40; 
                end
                
                C0HA2 = 5.84;
				C1HA2 = -0.02;
				C2HA2 = 5.42;
		  		dhA2 = -12.00;
							   
				C0HA3 = 6.05;
				C1HA3 = -0.07;
				C2HA3 = 5.55;
				dhA3 = 4.90;
					 		   
				C0C2 = 3.48; 
				C1C2 = -1.39; 
				C2C2 = 3.40; 
				dC2 = 2.80; 
					 		   
				C0C3 = 3.68; 
				C1C3 = -1.41; 
				C2C3 = 3.58; 
				dC3 = 15.60;
            elseif y(11) == 5
                if prompt_avg == 1
                    C0N = -3.25;
                    C1N = -0.20;
                    C2N = -3.19;
                    dN = -1.05;
				  
                elseif prompt_avg == 2
                    C0N = -2.56;
                    C1N = 0.38;
                    C2N = -2.52;
                    dN = -4.03;
                  
                else
                    C0N = -2.91;
                    C1N = 0.09;
                    C2N = -2.85;
                    dN = -2.54;
                  
                end
                C0HA = 6.82;
				C1HA = 0.22;
				C2HA = 6.22;
				dHA = -3.00;
                
                C0C = 4.12;
				C1C = -1.04;
				C2C = 4.00;
			 	dC = 7.08;
                
            elseif y(11) == 7
                if prompt_avg == 1
                    C0N = -2.20;
                    C1N = -0.55;
                    C2N = -2.73;
                    dN = -8.60;
				  
                    C0C = 3.48;
                    C1C = -1.39;
                    C2C = 3.40;
                	dC = 1.40;
                  
                elseif prompt_avg == 2
                    C0N = -1.74;
                    C1N = 0.20;
                    C2N = -1.86;
                    dN = -12.70;
				  
                    C0C = 3.18;
                    C1C = -1.69;
                    C2C = 3.08;
                    dC = -0.70;
                    
                else
                    C0N = -1.97;
                    C1N = -0.18;
                    C2N = -2.30;
                    dN = -10.65;
				  
                    C0C = 3.33;
                    C1C = -1.54;
                    C2C = 3.24;
                    dC = 0.35;
                    
                end
                C0HA = 5.50;
                C1HA = 0.53;
				C2HA = 4.99;
				dHA = -9.2;
                
            elseif y(11) == 3
                if prompt_avg == 1
                    C0N = -2.73;
                    C1N = -0.81;
                    C2N = -2.85;
                    dN = -1.2;
				  
                    C0C = 3.89;
                    C1C = -0.91;
                    C2C = 3.88;
                    dC = 8.2;
                
                elseif prompt_avg == 2
                    C0N = -2.25;
                    C1N = 0.09;
                    C2N = -2.22;
                    dN = -4.00;
				  
                    C0C = 3.76;
                    C1C = -1.30;
                    C2C = 3.74;
                    dC = 7.3;
                  
                else
                    C0N = -2.49;
                    C1N = -0.36;
                    C2N = -2.54;
                    dN = -2.6;
				  
                    C0C = 3.83;
                    C1C = -1.11;
                    C2C = 3.81;
                    dC = 7.75;
                  
                end
                C0HA = 6.33;
				C1HA = 0.68;
				C2HA = 5.69;
				dHA = -2.3;
                
            else
                if prompt_avg == 1
                    C0N2 = -3.25;
                    C1N2 = -0.20;
                    C2N2 = -3.19;
                    dN2 = -1.05;
                    
                elseif prompt_avg == 2
                    C0N2 = -2.56;
                    C1N2 = 0.38;
                    C2N2 = -2.52;
                    dN2 = -4.03;
                    
                elseif prompt_avg == 3
                    C0N2 = -2.91;
                    C1N2 = 0.09;
                    C2N2 = -2.85;
                    dN2 = -2.54;
                    
                end
                
                C0HA2 = 6.82;
                C1HA2 = 0.22;
                C2HA2 = 6.22;
                dHA2 = -3.00;
                
                C0C2 = 4.12;
                C1C2 = -1.04;
                C2C2 = 4.00;
                dC2 = 7.08;
                
                C0HA3 = C0HA2;
                C1HA3 = C1HA2;
                C2HA3 = C2HA2;
                dHA3 = dHA2;

                C0N3 = C0N2;
                C1N3 = C1N2;
                C2N3 = C2N2;
                dN3 = dN2;

                C0C3 = C0C2;
                C1C3 = C1C2;
                C2C3 = C2C2;
                dC3 = dC2;
            end
            
            %Phase Shift Section
            if y(13) == 1 
		 		dHA2 = 0;
		 		dN2 = 0; 
		 		dC2 = 0;
                
				dHA3 = 0; 
				dN3 = 0; 
				dC3 = 0;
                
                %In the case of beta branch
                dHA = 0;
		 		dN = 0; 
		 		dC = 0;
            end
            
            %Case variable setup
        end
        %beta branch calculation
        e = 0;
        clear d
        clear g
        if prompt_model == 3
            dHA_rad = dHA*pi/180;
            dN_rad = dN*pi/180;
            dC_rad = dC*pi/180;

            residue = y(1)
			for a=0:100
                for b=0:100-a
					c=100-b-a;
					e = e+1;
					
					% gauche negative
					GN_xrad = pi/3;
					GN_JHaHb = C0HA + C1HA*(cos((2*pi/3 - GN_xrad)-dHA_rad)) + C2HA*cos(2*((2*pi/3 - GN_xrad)-dHA_rad));
		    		GN_JNHb = C0N + C1N*(cos((-2*pi/3 - GN_xrad)-dN_rad)) + C2N*cos(2*((-2*pi/3 - GN_xrad)-dN_rad));
		 			GN_JCOHb = C0C + C1C*(cos((-GN_xrad)-dC_rad)) + C2C*cos(2*((-GN_xrad)-dC_rad));
					
					% trans 
					T_xrad = pi;
					T_JHaHb = C0HA + C1HA*(cos((2*pi/3 - T_xrad)-dHA_rad)) + C2HA*cos(2*((2*pi/3 - T_xrad)-dHA_rad));
		    		T_JNHb = C0N + C1N*(cos((-2*pi/3 - T_xrad)-dN_rad)) + C2N*cos(2*((-2*pi/3 - T_xrad)-dN_rad));
		 			T_JCOHb = C0C + C1C*(cos((-T_xrad)-dC_rad)) + C2C*cos(2*((-T_xrad)-dC_rad));

					% gauche pos  
					GP_xrad = -pi/3;
					GP_JHaHb = C0HA + C1HA*(cos((2*pi/3 - GP_xrad)-dHA_rad)) + C2HA*cos(2*((2*pi/3 - GP_xrad)-dHA_rad));
		    		GP_JNHb = C0N + C1N*(cos((-2*pi/3 - GP_xrad)-dN_rad)) + C2N*cos(2*((-2*pi/3 - GP_xrad)-dN_rad));
		 			GP_JCOHb = C0C + C1C*(cos((-GP_xrad)-dC_rad)) + C2C*cos(2*((-GP_xrad)-dC_rad));

					if y(11) == 3 
						% gauche neg 
				    	GN_JHaHb = C0HA + C1HA*(cos((-GN_xrad)-dHA_rad)) + C2HA*cos(2*((-GN_xrad)-dHA_rad));
				 		GN_JNHb = C0N + C1N*(cos((2*pi/3 - GN_xrad)-dN_rad)) + C2N*cos(2*((2*pi/3 - GN_xrad)-dN_rad));
						GN_JCOHb = C0C + C1C*(cos((-2*pi/3 - GN_xrad)-dC_rad)) + C2C*cos(2*((-2*pi/3 - GN_xrad)-dC_rad));
							
						% trans 	    		
						T_JHaHb=C0HA + C1HA*(cos((-T_xrad)-dHA_rad)) + C2HA*cos(2*((-T_xrad)-dHA_rad));
		 				T_JNHb=C0N + C1N*(cos((2*pi/3 - T_xrad)-dN_rad)) + C2N*cos(2*((2*pi/3 - T_xrad)-dN_rad));
						T_JCOHb=C0C + C1C*(cos((-2*pi/3 - T_xrad)-dC_rad)) + C2C*cos(2*((-2*pi/3 - T_xrad)-dC_rad));
							
						% gauche pos  
						GP_JHaHb=C0HA + C1HA*(cos((-GP_xrad)-dHA_rad)) + C2HA*cos(2*((-GP_xrad)-dHA_rad));
		 				GP_JNHb=C0N + C1N*(cos((2*pi/3 - GP_xrad)-dN_rad)) + C2N*cos(2*((2*pi/3 - GP_xrad)-dN_rad));
						GP_JCOHb=C0C + C1C*(cos((-2*pi/3 - GP_xrad)-dC_rad)) + C2C*cos(2*((-2*pi/3 - GP_xrad)-dC_rad));
					end 
					
					JHaHb = (a/100)*GN_JHaHb + (b/100)*T_JHaHb + (c/100)*GP_JHaHb;
					JNHb  = (a/100)*GN_JNHb  + (b/100)*T_JNHb  + (c/100)*GP_JNHb;
					JCOHb = (a/100)*GN_JCOHb + (b/100)*T_JCOHb + (c/100)*GP_JCOHb;
					
					TOCSYHbcalc = (sin(pi*JHaHb*TOCSYmix))^2;
			    	Int(1) = 1;
			    
			    	COHbcalc = (sin(pi*JCOHb*HNCOHBmix))^2;
			    	Int(2) = 1.836;

			    	NHbcalc = (sin(pi*JNHb*HNHBmix))^2;
			    	Int(3) = 0.6275;
			    	        
			    	Int(4) = (y(2)*Int(1)*TOCSYHbcalc + y(5)*Int(2)*COHbcalc + y(8)*Int(3)*NHbcalc) / (Int(1)^2*TOCSYHbcalc^2 + Int(2)^2*COHbcalc^2 + Int(3)^2* NHbcalc^2);
					
					Chi2 = ((y(2)) - (Int(4)*Int(1)*TOCSYHbcalc))^2/y(4) + ((y(5)) - (Int(4)*Int(2)*COHbcalc))^2/y(7) +  ((y(8)) - (Int(4)*Int(3)*NHbcalc))^2/y(10);
					% the constants above were estimated from model 1 residues, corresponding to the reciprocal of the max intensities observable for the three spectra, HNHB was the most intense, HNCOHB was the weakest 
					
					min1 = Chi2;
		    			
		    		g(e,1) = y(1);
					g(e,2) = min1;
					g(e,3) = Int(4)*Int(1)*TOCSYHbcalc;
					g(e,4) = Int(4)*Int(2)*COHbcalc;
					g(e,5) = Int(4)*Int(3)*NHbcalc;
					g(e,6) = a;
					g(e,7) = b;
					g(e,8) = c;
				
					a1 = a + 1;
					b1 = b + 1;
					L(a1,b1) = min1; 
		
                end
            end	  
        %Create Table for Betabranch    
		[M,I] = min(g);
		final_int = I(2);
        
        TBetaBranch(counter + 1, :) = g(final_int, :);
        
        for xTest=0:101
            for yTest=0:101
                if xTest + yTest >= 102
                   L(xTest,yTest) = NaN; 
                end
            end
        end
		%contourf(t,z,L,50); % create contour plot
		[c,H] = contourf(t,z,L,20);
        
        box off
        set(gca, 'TickDir','out');
        
        fontname(fontType)

        colourbarTitle = contourcbar('westoutside');
        title(colourbarTitle,'Χ^{2}','FontSize',fontSize)
        
        graphTitle = sprintf('Beta Branch Contour Plot for Residue %d',y(1));
        title(graphTitle,'FontSize',fontSize)
        xlabel('% Gauche Negative','FontSize',fontSize)
        ylabel('% Trans','FontSize',fontSize)
        
        fig = gcf;
        
        betaBranchFileName = sprintf('BetaBranch%d',y(1));
        betaBranchFull = fullfile(betaBranchDir,betaBranchFileName);
        betaBranchFullName = sprintf('%s.png',betaBranchFull);
        exportgraphics(fig,betaBranchFullName,'Resolution',300)
        
        %Flexible Calculation
        else
            dHA_rad2 = dHA2*pi/180;
            dN_rad2 = dN2*pi/180;
            dC_rad2 = dC2*pi/180;
            dHA_rad3 = dHA3*pi/180;
            dN_rad3 = dN3*pi/180;
            dC_rad3 = dC3*pi/180;
              
            residue = y(1)
            for a=0:100  %a is the proportion of gauche negative
                for b=0:100-a    %b is the proportion of trans
				    c=100-b-a;	%c is the proportion of gauche positive
				    e = e+1;	%e is the counter for the table with all the values of a,b,c

				    % gauche negative
				    GN_xrad = pi/3;
				    GN_JHaHb2 = C0HA2 + C1HA2*(cos((2*pi/3 - GN_xrad)-dHA_rad2)) + C2HA2*cos(2*((2*pi/3 - GN_xrad)-dHA_rad2));
				    GN_JHaHb3 = C0HA3 + C1HA3*(cos((-GN_xrad)-dHA_rad3)) + C2HA3*cos(2*((-GN_xrad)-dHA_rad3));
				    GN_JNHb2 = C0N2 + C1N2*(cos((-2*pi/3 - GN_xrad)-dN_rad2)) + C2N2*cos(2*((-2*pi/3 - GN_xrad)-dN_rad2));
				    GN_JNHb3 = C0N3 + C1N3*(cos((2*pi/3 - GN_xrad)-dN_rad3)) + C2N3*cos(2*((2*pi/3 - GN_xrad)-dN_rad3));
				    GN_JCOHb2 = C0C2 + C1C2*(cos((-GN_xrad)-dC_rad2)) + C2C2*cos(2*((-GN_xrad)-dC_rad2));
				    GN_JCOHb3 = C0C3 + C1C3*(cos((-2*pi/3 - GN_xrad)-dC_rad3)) + C2C3*cos(2*((-2*pi/3 - GN_xrad)-dC_rad3));
						
				    % trans 
				    T_xrad = pi;
				    T_JHaHb2 = C0HA2 + C1HA2*(cos((2*pi/3 - T_xrad)-dHA_rad2)) + C2HA2*cos(2*((2*pi/3 - T_xrad)-dHA_rad2));
				    T_JHaHb3 = C0HA3 + C1HA3*(cos((-T_xrad)-dHA_rad3)) + C2HA3*cos(2*((-T_xrad)-dHA_rad3));
				    T_JNHb2 = C0N2 + C1N2*(cos((-2*pi/3 - T_xrad)-dN_rad2)) + C2N2*cos(2*((-2*pi/3 - T_xrad)-dN_rad2));
				    T_JNHb3 = C0N3 + C1N3*(cos((2*pi/3 - T_xrad)-dN_rad3)) + C2N3*cos(2*((2*pi/3 - T_xrad)-dN_rad3));
				    T_JCOHb2 = C0C2 + C1C2*(cos((-T_xrad)-dC_rad2)) + C2C2*cos(2*((-T_xrad)-dC_rad2));
				    T_JCOHb3 = C0C3 + C1C3*(cos((-2*pi/3 - T_xrad)-dC_rad3)) + C2C3*cos(2*((-2*pi/3 - T_xrad)-dC_rad3));
						
				    % gauche pos  
				    GP_xrad = -pi/3;
				    GP_JHaHb2 = C0HA2 + C1HA2*(cos((2*pi/3 - GP_xrad)-dHA_rad2)) + C2HA2*cos(2*((2*pi/3 - GP_xrad)-dHA_rad2));
			    	GP_JHaHb3 = C0HA3 + C1HA3*(cos((-GP_xrad)-dHA_rad3)) + C2HA3*cos(2*((-GP_xrad)-dHA_rad3));
			    	GP_JNHb2 = C0N2 + C1N2*(cos((-2*pi/3 - GP_xrad)-dN_rad2)) + C2N2*cos(2*((-2*pi/3 - GP_xrad)-dN_rad2));
			 	    GP_JNHb3 = C0N3 + C1N3*(cos((2*pi/3 - GP_xrad)-dN_rad3)) + C2N3*cos(2*((2*pi/3 - GP_xrad)-dN_rad3));
			 	    GP_JCOHb2 = C0C2 + C1C2*(cos((-GP_xrad)-dC_rad2)) + C2C2*cos(2*((-GP_xrad)-dC_rad2));
				    GP_JCOHb3 = C0C3 + C1C3*(cos((-2*pi/3 - GP_xrad)-dC_rad3)) + C2C3*cos(2*((-2*pi/3 - GP_xrad)-dC_rad3));
										
				    JHaHb2 = (a/100)*GN_JHaHb2 + (b/100)*T_JHaHb2 + (c/100)*GP_JHaHb2;
				    JHaHb3 = (a/100)*GN_JHaHb3 + (b/100)*T_JHaHb3 + (c/100)*GP_JHaHb3;
				    JNHb2  = (a/100)*GN_JNHb2  + (b/100)*T_JNHb2  + (c/100)*GP_JNHb2;
				    JNHb3  = (a/100)*GN_JNHb3  + (b/100)*T_JNHb3  + (c/100)*GP_JNHb3;
				    JCOHb2 = (a/100)*GN_JCOHb2 + (b/100)*T_JCOHb2 + (c/100)*GP_JCOHb2;
				    JCOHb3 = (a/100)*GN_JCOHb3 + (b/100)*T_JCOHb3 + (c/100)*GP_JCOHb3;
					
				    TOCSYHb2calc = (sin(pi*JHaHb2*TOCSYmix))^2;
			    	TOCSYHb3calc = (sin(pi*JHaHb3*TOCSYmix))^2;

			    	Int(1) = (TOCSYHb2calc*y(2) + TOCSYHb3calc*y(3)) / (TOCSYHb2calc^2 + TOCSYHb3calc^2);
			    
			    	COHb2calc = (sin(pi*JCOHb2*HNCOHBmix))^2;
			    	COHb3calc = (sin(pi*JCOHb3*HNCOHBmix))^2;
			    
                    Int(2) = (COHb2calc*y(5) + COHb3calc*y(6)) / (COHb2calc^2 + COHb3calc^2);

                    NHb2calc = (sin(pi*JNHb2*HNHBmix))^2;
                    NHb3calc = (sin(pi*JNHb3*HNHBmix))^2;

                    Int(3) = (NHb2calc*y(8) + NHb3calc*y(9)) / (NHb2calc^2 + NHb3calc^2);

                    Chi2 = ((y(2) - Int(1)*TOCSYHb2calc)^2)/y(4) + ((y(3) - Int(1)*TOCSYHb3calc)^2)/y(4) + ((y(5) - Int(2)*COHb2calc)^2)/y(7) + ((y(6) - Int(2)*COHb3calc)^2)/y(7) +  ((y(8) - Int(3)*NHb2calc)^2)/y(10) + ((y(9) - Int(3)*NHb3calc)^2)/y(10);

                    Int(4) = 1;
				        		    
				    min1 = Chi2; 

					if prompt_model == 4;
                        Int(1) = 1; %For model 4, degenerate HB2/HB3, set Int(1), Int(2), and Int(3) as
						Int(2) = 1.836;% relative intensities for HTOCSY, HN(CO)HB, and HNHB spectra
						Int(3) = 0.6275; % then Int(4) becomes the scaling factor for all 3 spectra
						Int(4) = (y(2)*Int(1)*(TOCSYHb2calc + TOCSYHb3calc) + y(5)*Int(2)*(COHb2calc + COHb3calc) + y(8)*Int(3)*(NHb2calc + NHb3calc)) / (Int(1)^2*(TOCSYHb2calc+TOCSYHb3calc)^2 + Int(2)^2*(COHb2calc+COHb3calc)^2+ Int(3)^2*(NHb2calc+NHb3calc)^2);
						
						Chi2 = (y(2) - (Int(1)*Int(4)*(TOCSYHb2calc+TOCSYHb3calc)))^2/y(4) + (y(5) - (Int(2)*Int(4)*(COHb2calc+COHb3calc)))^2/y(7) +  (y(8) - (Int(3)*Int(4)*(NHb2calc+NHb3calc)))^2/y(10); 
						
						min1 = Chi2;
					end 
					d(e,1) = y(1);
					d(e,2) = min1;
					d(e,3) = TOCSYHb2calc*Int(1)*Int(4);
					d(e,4) = TOCSYHb3calc*Int(1)*Int(4);
					d(e,5) = COHb2calc*Int(2)*Int(4);
					d(e,6) = COHb3calc*Int(2)*Int(4);
					d(e,7) = NHb2calc*Int(3)*Int(4);
					d(e,8) = NHb3calc*Int(3)*Int(4);
					d(e,9) = a;
					d(e,10) = b;
					d(e,11) = c;
					
					a1 = a + 1;
					b1 = b + 1;
					L(a1,b1) = min1; 
						
                end
            end	
			
            %Create table for flexible/degenerate
            [M,I] = min(d);
            final_int = I(2);
            
            TFlexOrDegen(counter + 1,:) = d(final_int,:);
            
            
            for xTest=0:101
                for yTest=0:101
                    if xTest + yTest >= 102
                       L(xTest,yTest) = NaN; 
                    end
                end
            end
            %contourPlot = contourf(t,z,L,50); % create contour plot 
            [c,H] = contourf(t,z,L,20);
            
            box off
            set(gca, 'TickDir','out');

            fontname(fontType)
            
            colourbarTitle = contourcbar('westoutside');
            title(colourbarTitle,'Χ^{2}','FontSize',fontSize)
            
            graphTitle = sprintf('Flexible/Degnerate Contour Plot for Residue %d',y(1));
            title(graphTitle,'FontSize',fontSize)
            xlabel('% Gauche Negative','FontSize',fontSize)
            ylabel('% Trans','FontSize',fontSize)
            
            fig = gcf;
        
            flexibleDegenFileName = sprintf('FlexibleDegen%d', y(1));
            flexibleDegenFull = fullfile(flexibleDegenDir,flexibleDegenFileName);
            flexibleDegenFullName = sprintf('%s.png',flexibleDegenFull);
            exportgraphics(fig,flexibleDegenFullName,'Resolution',300)

            %Runs flexible residues through the rigid calculation
            if prompt_model == 1;
			% The variable x is the chi1 dihedral angle, defined according to IUPAC as looking down the Ca-Cb bond and measuring the dihedral angle between N and Cg, with a positive dihedral being a clockwise rotation of Cg relative to N
                for x=1:360
                    xrad = x*pi/180;
                    JHaHb2=C0HA2 + C1HA2*(cos((2*pi/3 - xrad)-dHA_rad2)) + C2HA2*cos(2*((2*pi/3 - xrad)-dHA_rad2));
                    JHaHb3=C0HA3 + C1HA3*(cos((-xrad)-dHA_rad3)) + C2HA3*cos(2*((-xrad)-dHA_rad3));
                    JNHb2=C0N2 + C1N2*(cos((-2*pi/3 - xrad)-dN_rad2)) + C2N2*cos(2*((-2*pi/3 - xrad)-dN_rad2));
                    JNHb3=C0N3 + C1N3*(cos((2*pi/3 - xrad)-dN_rad3)) + C2N3*cos(2*((2*pi/3 - xrad)-dN_rad3));
                    JCOHb2=C0C2 + C1C2*(cos((-xrad)-dC_rad2)) + C2C2*cos(2*((-xrad)-dC_rad2));
                    JCOHb3=C0C3 + C1C3*(cos((-2*pi/3 - xrad)-dC_rad3)) + C2C3*cos(2*((-2*pi/3 - xrad)-dC_rad3));			     

                    TOCSYHb2calc = (sin(pi*JHaHb2*TOCSYmix))^2;
                    TOCSYHb3calc = (sin(pi*JHaHb3*TOCSYmix))^2;

                    Int(1) = (TOCSYHb2calc*y(2) + TOCSYHb3calc*y(3)) / (TOCSYHb2calc^2 + TOCSYHb3calc^2);

                    COHb2calc = (sin(pi*JCOHb2*HNCOHBmix))^2;
                    COHb3calc = (sin(pi*JCOHb3*HNCOHBmix))^2;

                    Int(2) = (COHb2calc*y(5) + COHb3calc*y(6)) / (COHb2calc^2 + COHb3calc^2);

                    NHb2calc = (sin(pi*JNHb2*HNHBmix))^2;
                    NHb3calc = (sin(pi*JNHb3*HNHBmix))^2;

                    Int(3) = (NHb2calc*y(8) + NHb3calc*y(9)) / (NHb2calc^2 + NHb3calc^2);

                    Chi2 = ((y(2) - Int(1)*TOCSYHb2calc)^2)/y(4) + ((y(3) - Int(1)*TOCSYHb3calc)^2)/y(4) + ((y(5) - Int(2)*COHb2calc)^2)/y(7) + ((y(6) - Int(2)*COHb3calc)^2)/y(7) +  ((y(8) - Int(3)*NHb2calc)^2)/y(10) + ((y(9) - Int(3)*NHb3calc)^2)/y(10);
                    
                    if ( x >= 0 && x <= 30 ) || ( x >= 90 && x <= 150 ) || ( x >= 210 && x <= 270 ) || ( x >= 330 && x <= 360 )
                        value = Chi2 + 10;
                    else    
                        value = Chi2;
                    end

                    a3(x,1) = residue;
                    a3(x,2) = x;
                    a3(x,3) = value;
                    a3(x,4) = Int(1);% theoretical max intensity for HTOCSY in this residue
                    a3(x,5) = Int(2);% theoretical max intensity for HN(CO)HB in this residue
                    a3(x,6) = Int(3);% theoretical max intensity for HNHB in this residue
                    a3(x,7) = Int(1)*TOCSYHb2calc;
                    a3(x,8) = Int(1)*TOCSYHb3calc;
                    a3(x,9) = Int(2)*COHb2calc;
                    a3(x,10) = Int(2)*COHb3calc;
                    a3(x,11) = Int(3)*NHb2calc; %#ok<SAGROW>
                    a3(x,12) = Int(3)*NHb3calc;

                end

                %Create table for rigid
                [M,I] = min(a3); %This generates a vector M with minima for each column of a, vector I is the index or each minimum value
                finalchi1degrees = I(3); %finalchi1 is the chi1 dihedral angle in degrees that has the minimum chi2 value
                
                TRigid(counter + 1,:) = a3(finalchi1degrees,:);

                plot(a3(:,2),a3(:,3))
                xlim([0,360])
                set(gca,'xtick',[0:60:360]);
                set(gca, 'TickDir','in');

                fontname(fontType)
                
                graphTitle = sprintf('Rigid Structure Prediction for Residue %d',y(1));
                title(graphTitle,'FontSize',fontSize)
                xlabel('Χ_{1} Angle (°)','FontSize',fontSize)
                ylabel('Χ^{2}','FontSize',fontSize)
                
                
                fig = gcf;
        
                rigidFileName = sprintf('Rigid%d',y(1));
                rigidFull = fullfile(rigidDir,rigidFileName);
                rigidFullName = sprintf('%s.png',rigidFull);
                exportgraphics(fig,rigidFullName,'Resolution',300)

                yold(2)=y(2);
                yold(3)=y(3);
                yold(5)=y(5);
                yold(6)=y(6);
                yold(8)=y(8);
                yold(9)=y(9);

                for errorloop = 1:20;
                    y(2) = yold(2) + normrnd(0,1)*y(4);
                    y(3) = yold(3) + normrnd(0,1)*y(4);
                    y(5) = yold(5) + normrnd(0,1)*y(7);
                    y(6) = yold(6) + normrnd(0,1)*y(7);
                    y(8) = yold(8) + normrnd(0,1)*y(10);
                    y(9) = yold(9) + normrnd(0,1)*y(10);

                    for x=1:360
                        xrad = x*pi/180;
                        JHaHb2=C0HA2 + C1HA2*(cos((2*pi/3 - xrad)-dHA_rad2)) + C2HA2*cos(2*((2*pi/3 - xrad)-dHA_rad2));
                        JHaHb3=C0HA3 + C1HA3*(cos((-xrad)-dHA_rad3)) + C2HA3*cos(2*((-xrad)-dHA_rad3));
                        JNHb2=C0N2 + C1N2*(cos((-2*pi/3 - xrad)-dN_rad2)) + C2N2*cos(2*((-2*pi/3 - xrad)-dN_rad2));
			            JNHb3=C0N3 + C1N3*(cos((2*pi/3 - xrad)-dN_rad3)) + C2N3*cos(2*((2*pi/3 - xrad)-dN_rad3));
                        JCOHb2=C0C2 + C1C2*(cos((-xrad)-dC_rad2)) + C2C2*cos(2*((-xrad)-dC_rad2));
                        JCOHb3=C0C3 + C1C3*(cos((-2*pi/3 - xrad)-dC_rad3)) + C2C3*cos(2*((-2*pi/3 - xrad)-dC_rad3));

                        TOCSYHb2calc = (sin(pi*JHaHb2*TOCSYmix))^2;
                        TOCSYHb3calc = (sin(pi*JHaHb3*TOCSYmix))^2;

                        Int(1) = (TOCSYHb2calc*y(2) + TOCSYHb3calc*y(3)) / (TOCSYHb2calc^2 + TOCSYHb3calc^2);

                        COHb2calc = (sin(pi*JCOHb2*HNCOHBmix))^2;
                        COHb3calc = (sin(pi*JCOHb3*HNCOHBmix))^2;

                        Int(2) = (COHb2calc*y(5) + COHb3calc*y(6)) / (COHb2calc^2 + COHb3calc^2);

                        NHb2calc = (sin(pi*JNHb2*HNHBmix))^2;
                        NHb3calc = (sin(pi*JNHb3*HNHBmix))^2;

                        Int(3) = (NHb2calc*y(8) + NHb3calc*y(9)) / (NHb2calc^2 + NHb3calc^2);

                        Chi2 = ((y(2) - Int(1)*TOCSYHb2calc)^2)/y(4) + ((y(3) - Int(1)*TOCSYHb3calc)^2)/y(4) + ((y(5) - Int(2)*COHb2calc)^2)/y(7) + ((y(6) - Int(2)*COHb3calc)^2)/y(7) +  ((y(8) - Int(3)*NHb2calc)^2)/y(10) + ((y(9) - Int(3)*NHb3calc)^2)/y(10);

                        value = Chi2; 

                        if x >= 0 && x <= 20; 
                            value = Chi2 + 10;
                        end 

                        if x >= 100 && x <= 140; 
                            value = Chi2 + 10;
                        end 

                        if x >= 220 && x <= 260; 
                            value = Chi2 + 10;
                        end 

                        if x >= 340 && x <= 360; 
                            value = Chi2 + 10;
                        end 
                        
                        a3(x,1) = residue;
                        a3(x,2) = x;
                        a3(x,3) = value;
                        a3(x,4) = Int(1); % theoretical max intensity for HTOCSY in this residue
                        a3(x,5) = Int(2);
                        a3(x,6) = Int(3);
                        a3(x,7) = Int(1)*TOCSYHb2calc;
                        a3(x,8) = Int(1)*TOCSYHb3calc;
                        a3(x,9) = Int(2)*COHb2calc;
                        a3(x,10) = Int(2)*COHb3calc;
                        a3(x,11) = Int(3)*NHb2calc;
                        a3(x,12) = Int(3)*NHb3calc;

                    end

                    [M,I] = min(a3); %This generates a vector M with minima for each column of a, vector I is the index or each minimum value
                    finalchi1degrees = I(3); %finalchi1 is the chi1 dihedral angle in degrees that has the minimum chi2 value
                    
                    % finalchi1= finalchi1degrees*pi/180;
                    %finaltab has the following columns: residue chi1 chi2 JHaHb2 JHaHb3 JNHb2 JNHb3 JCOHb2 JCOHb3 TOCSYmaxintensity HNCOHBmaxintensity HNHBmaxintensity
                    MC_chi1(errorloop)=a3(finalchi1degrees,2);

                end
                MC_std = std(MC_chi1);
                TRigid(counter + 1,13) = MC_std;
            end
        end
    end
    counter = counter + 1;
end

writematrix(TBetaBranch,'BetaBranch.txt','Delimiter',' ')
writematrix(TFlexOrDegen,'FlexOrDegen.txt','Delimiter',' ')
writematrix(TRigid,'Rigid.txt','Delimiter',' ')