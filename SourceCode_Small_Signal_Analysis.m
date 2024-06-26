%1. Check the format and units of TXT dataset.
%	Time(s) \t Voltage(V) \t Current(mA)
	% For BioLogic instruments, the default current unit is mA. 
		% it will be converted to A unit. 
	% otherwise, edit the code in 'Convert Unit (mA --> A)' section 

%2. Edit the file name; 
	% "Raw_GS.txt" and "Raw_DS.txt" for the dataset from gate and drain electrodes.

%3. Assign input parameters; 
	clear Para Dim
	% Assign dimension of active layer: 
	Dim.L_ch = (20E-6)*100; 		% [cm] Channel length
	Dim.w = (100E-6)*100; 			% [cm] Channel Width 
	Dim.d = (36.12E-9)*100; 			% [cm] Channel thickness
	Dim.L_contact = (100E-6)*100; 	% [cm] Contact length (we considered 'symmetric structure' of active layer.)
	Dim.A_active = (Dim.L_ch+2*Dim.L_contact)*(Dim.w); 
									% [cm^3] volume of active material
	Dim.A_ch = Dim.L_ch*Dim.w; 
									% [cm^3] volume of channel area    
	Dim.ArealFactor = Dim.A_ch/Dim.A_active; 
									% Area correction factor 
	% Assign measurement conditions: 
	Para.dTs = 0.001; 				% [sec] Sampling Rate
	Para.f_AC = 10; 				% [Hz] Input gate voltage frequency
	Para.Freq_Filter = [0.1 2]; 	% [#1, #2] BandPassFilter; Frequency Ratio (eg. #1*Freq ~ #2*Freq)
	Para.FFTCycle = 1; 				% [#]  Number of cycles for Fourier transform calculation
	Para.Fs = 1/(Para.dTs);			% [Hz] Sampling Frequency


%4. Import raw data (TXT format)
	warning off
	opts = delimitedTextImportOptions("NumVariables", 3);
	opts.DataLines = [2, Inf]; 
	opts.Delimiter = "\t";
	opts.VariableNames = ["time[s]", "V[V]", "I[mA]"];
	opts.VariableTypes = ["double", "double", "double"];
	opts.ExtraColumnsRule = "ignore";
	opts.EmptyLineRule = "read";
	Raw_GS = rmmissing(table2array(readtable('Raw_GS', opts))); % Remove empty rows 
	Raw_DS = rmmissing(table2array(readtable('Raw_DS', opts))); % Remove empty rows 
	% Convert Unit (mA --> A)
		Raw_GS(:,3) = Raw_GS(:,3)*1E-3;
		Raw_DS(:,3) = Raw_DS(:,3)*1E-3;
	clear opts
	warning on
	% Assign Region of Interest
				%% below is to assign whole range of dataset. 
				%% to assign specific range, edit the 'Para.InitPoint' and 'Para.EndPoint' . 
	Time_Margine = 10000; 			% Time margin for better signal filtering. 
	Para.InitPoint = Time_Margine+1;
									% [Row #] of the begining of ROI
	Para.EndPoint = min(length(Raw_GS), length(Raw_DS))-Time_Margine;
									% [Row #] of the end of ROI

%5. Assign the ROI.
	Para.InitPoint = Para.InitPoint-Time_Margine;
	Para.EndPoint = Para.EndPoint+Time_Margine;
	Temp.Data_GS = Raw_GS(Para.InitPoint:Para.EndPoint,:);
	Temp.Data_DS = Raw_DS(Para.InitPoint:Para.EndPoint,:);
	% assign 'time' from 'row number'. 
	RoiTimeFrom = max((Temp.Data_GS(1,1)),(Temp.Data_DS(1,1)));
	RoiTimeEnd = min((Temp.Data_GS(end,1)),(Temp.Data_DS(end,1)));
	RoiTime = [RoiTimeFrom RoiTimeEnd]; % [s]
	clear RoiTimeFrom RoiTimeEnd

%6. interpolation of dataset (to assign uniform time interval)
	Temp.Data_GS = unique(rmmissing(Temp.Data_GS),'rows','sorted'); %sort as time
	A_ITP_GS(:,1) = RoiTime(1):Para.dTs:RoiTime(2); %define new time scale
	for i = 2:width(Temp.Data_GS) %interpolation for all row based on time.
		A_ITP_GS(:,i) = interp1(Temp.Data_GS(:,1),Temp.Data_GS(:,i),A_ITP_GS(:,1));     
	end
	A_ITP_GS(isnan(A_ITP_GS)) = 0; %convert 'NaN' to '0'
	Temp.Data_DS = unique(rmmissing(Temp.Data_DS),'rows','sorted'); %sort as time
	A_ITP_DS(:,1) = RoiTime(1):Para.dTs:RoiTime(2); %define new time scale
	for i = 2:width(Temp.Data_DS) %interpolation for all row based on time.
		A_ITP_DS(:,i) = interp1(Temp.Data_DS(:,1),Temp.Data_DS(:,i),A_ITP_GS(:,1));     
	end
	A_ITP_DS(isnan(A_ITP_DS)) = 0; %convert 'NaN' to '0'

%7. Assign ROI data. 
	VGS.Roi = A_ITP_GS(:,2);
	IGS.Roi = A_ITP_GS(:,3);
	VDS.Roi = A_ITP_DS(:,2);
	IDS.Roi = A_ITP_DS(:,3);
	Para.Time = A_ITP_GS(:,1);
    clear A_ITP_GS A_ITP_DS

%8. Remove high frequency noise with LPF
	CutOffFreq = Para.f_AC*Para.Freq_Filter(2); 
	VGS.Roi = VGS.Roi-highpass(VGS.Roi, CutOffFreq,Para.Fs);
	IGS.Roi = IGS.Roi-highpass(IGS.Roi, CutOffFreq,Para.Fs);
	VDS.Roi = VDS.Roi-highpass(VDS.Roi, CutOffFreq,Para.Fs);
	IDS.Roi = IDS.Roi-highpass(IDS.Roi, CutOffFreq,Para.Fs);
	clear CutOffFreq 

%9. Obtain AC component
	%Cut-off frequency = [input frequency] x [1/10]
	CutOffFreq = Para.f_AC*Para.Freq_Filter(1); 
	VGS.AC = highpass(VGS.Roi,CutOffFreq,Para.Fs);
	IGS.AC = highpass(IGS.Roi,CutOffFreq,Para.Fs);
	VDS.AC = highpass(VDS.Roi,CutOffFreq,Para.Fs);
	IDS.AC = highpass(IDS.Roi,CutOffFreq,Para.Fs);
	clear CutOffFreq 

%10. Obtain DC component
	VDS.DC = VDS.Roi-VDS.AC;
	IDS.DC = IDS.Roi-IDS.AC;
	VGS.DC = VGS.Roi-VGS.AC;
	IGS.DC = IGS.Roi-IGS.AC;
	Para.VGS = VGS.DC; 

%11. Assign parameter array with NAN;
	VGS.Ang(1:length(Para.Time),1) = nan;
	VDS.Ang(1:length(Para.Time),1) = nan;
	IGS.Ang(1:length(Para.Time),1) = nan;
	IDS.Ang(1:length(Para.Time),1) = nan;
	VGS.AC_A(1:length(Para.Time),1) = nan;
	VDS.AC_A(1:length(Para.Time),1) = nan;
	IGS.AC_A(1:length(Para.Time),1) = nan;
	IDS.AC_A(1:length(Para.Time),1) = nan;
	IGS.Phase(1:length(Para.Time),1) = nan;
	IDS.Phase(1:length(Para.Time),1) = nan;

%12. Fourier transform 
	Para.FT_Num = 1/(Para.dTs*Para.f_AC);
	npts = round(Para.FFTCycle*Para.FT_Num);
	nptsHalf = fix(npts/2);
	for i = 1:length(Para.Time)
		if i > (nptsHalf) && i<length(VGS.AC)-nptsHalf+1
			VGS.Token = fft(VGS.AC(i-nptsHalf:i+nptsHalf));
			IGS.Token = fft(IGS.AC(i-nptsHalf:i+nptsHalf));
			IDS.Token = fft(IDS.AC(i-nptsHalf:i+nptsHalf));
			[~, idx_VGS] = max(abs(VGS.Token));
			VGS.Ang(i,1) = (angle(VGS.Token(idx_VGS)));
			IGS.Ang(i,1) = (angle(IGS.Token(idx_VGS)));
			IDS.Ang(i,1) = (angle(IDS.Token(idx_VGS)));
			VGS.AC_A(i,1) = 2*abs(VGS.Token(idx_VGS))/npts;
			IGS.AC_A(i,1) = 2*abs(IGS.Token(idx_VGS))/npts;
			IDS.AC_A(i,1) = 2*abs(IDS.Token(idx_VGS))/npts;
			IGS.Phase(i,1) = (IGS.Ang(i) - VGS.Ang(i));
			IDS.Phase(i,1) = (IDS.Ang(i) - VGS.Ang(i));
		end  
	end
	clear idx_VGS npts nptsHalf NumUniquePts VGS.Token IGS.Token IDS.Token

%13. Define real and imaginary part from AC data. 
	IGS.AC_R = (IGS.AC_A.*cos(IGS.Phase));
	IGS.AC_I = (IGS.AC_A.*sin(IGS.Phase));
	% Remove a crosstalk between gate and drain current; Assume that the active layer is symmetric on the channel region.
	IDS.AC_R = (IDS.AC_A.*cos(IDS.Phase))+(IGS.AC_R/2);
	IDS.AC_I = (IDS.AC_A.*sin(IDS.Phase))+(IGS.AC_I/2);
	IDS.AC_A = ((IDS.AC_I.^2)+(IDS.AC_R.^2)).^0.5;

%14. Calculate the parameters
	% [C] Capacitance [F]
	Para.C.AC = (IGS.AC_I/(2*pi*Para.f_AC))./((VGS.AC_A));
	% [CVol] Volumetric capacitance [F cm-3]
	Para.CVol.AC = Para.C.AC./(Dim.A_active*Dim.d); 
	% [Gm] Transconductance = dIDS / dVGS [S]
	Para.Gm.AC = -IDS.AC_R/nanmedian(VGS.AC_A);
	Para.Gm.DC = -[(diff(IDS.DC)./diff(VGS.DC));0];
	% [G] Conductance = IDS/VDS [S]
	Para.G = IDS.DC/mean(VDS.DC); 
	% [Sigma] Conductivity [S/cm]; dimension-normalized conductance G
	Para.Sigma = Para.G/((Dim.w*Dim.d)/Dim.L_ch);
	% [Tau] TimeConstant [s]
	Para.Tau.AC = ((IGS.AC_I/(2*pi*Para.f_AC))*Dim.ArealFactor)./IDS.AC_R;
	% [Mob] Mobility [cm2 V-1s-1]
	Para.Mob.AC = -(Dim.L_ch^2)./(mean(VDS.DC)*Para.Tau.AC);
	% [uC] Mobility x volumetric capacitance [F cm-1 V-1 s-1]
	Para.uC.DC = Para.Gm.DC./(((Dim.w*Dim.d)/Dim.L_ch)*mean(VDS.DC));
	Para.uC.AC = Para.Mob.AC.*Para.CVol.AC;
	Para.Mob.DC = Para.uC.DC./Para.CVol.AC;




% Export data
% Assign row # to export. (from "ROI_Init" to "ROI_End")

ROI_Init=40000;
ROI_End=740000;
% File name; Result.txt; will be saved in identical folder. 
filenameT = sprintf('Result.txt');

	Data=[];
		Data(:,1)=VGS.DC(ROI_Init:ROI_End);
		Data(:,2)=VGS.Roi(ROI_Init:ROI_End);
		Data(:,3)=IGS.Roi(ROI_Init:ROI_End);
		Data(:,4)=VDS.Roi(ROI_Init:ROI_End);
		Data(:,5)=IDS.Roi(ROI_Init:ROI_End);
		Data(:,6)=VGS.DC(ROI_Init:ROI_End);
		Data(:,7)=IGS.DC(ROI_Init:ROI_End);
		Data(:,8)=VDS.DC(ROI_Init:ROI_End);
		Data(:,9)=IDS.DC(ROI_Init:ROI_End);
		Data(:,10)=VGS.AC(ROI_Init:ROI_End);
		Data(:,11)=IGS.AC(ROI_Init:ROI_End);
		Data(:,12)=IGS.AC_R(ROI_Init:ROI_End);
		Data(:,13)=IGS.AC_I(ROI_Init:ROI_End);
		Data(:,14)=VDS.AC(ROI_Init:ROI_End);
		Data(:,15)=IDS.AC(ROI_Init:ROI_End);
		Data(:,16)=IDS.AC_R(ROI_Init:ROI_End);
		Data(:,17)=IDS.AC_I(ROI_Init:ROI_End);
		Data(:,18)=Para.G(ROI_Init:ROI_End);
		Data(:,19)=Para.Sigma(ROI_Init:ROI_End);
		Data(:,20)=Para.C.AC(ROI_Init:ROI_End);
		Data(:,21)=Para.CVol.AC(ROI_Init:ROI_End);
		Data(:,22)=Para.Gm.AC(ROI_Init:ROI_End);
		Data(:,23)=Para.Gm.DC(ROI_Init:ROI_End);
		Data(:,24)=Para.uC.AC(ROI_Init:ROI_End);
		Data(:,25)=Para.uC.DC(ROI_Init:ROI_End);
		Data(:,26)=Para.Tau.AC(ROI_Init:ROI_End);
		Data(:,27)=Para.Mob.AC(ROI_Init:ROI_End);
		Data(:,28)=Para.Mob.DC(ROI_Init:ROI_End);
		Data(:,29)=Para.Time(ROI_Init:ROI_End);

	T = array2table(Data);
	TSize=size(T);
	T.Properties.VariableNames(1:TSize(2)) = {'VGS.DC[V]','VGS.Roi[V]','IGS.Roi[A]','VDS.Roi[V]','IDS.Roi[A]',...
											'VGS.DC_[V]','IGS.DC[A]','VDS.DC[V]','IDS.DC[A]'... 	
											'VGS.AC[V]','IGS.AC[A]','IGS.AC_R[A]','IGS.AC_I[A]',...
											'VDS.AC[V]','IDS.AC[A]','IDS.AC_R[A]','IDS.AC_I[A]'...
											'Para.G[S]','Para.Sigma[Scm-1]','Para.C.AC[F]','Para.CVol.AC[Fcm-3]',...
											'Para.Gm.AC[S]','Para.Gm.DC[S]','Para.uC.AC[Fcm-1V-1s-1]','Para.uC.DC[Fcm-1V-1s-1]',...
											'Para.Tau.AC[s]','Para.Mob.AC[cm2V-1s-1]','Para.Mob.DC[cm2V-1s-1]','Time[s]'};

	writetable(T,filenameT,'delimiter','tab');



	figure
	plot(Data(:,1),Data(:,27))
	hold on
	plot(Data(:,1),Data(:,28))
	title('mobility \mu')
	xlabel('V_{GS.DC} (V)')
	ylabel('\mu (cm^2 V^{-1} s^{-1})')
	legend({'\mu_{AC}','\mu_{DC}'})
	ylim([-0.5 3])
	hold off
