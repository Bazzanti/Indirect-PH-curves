%% -------------------------------------------------------------- preamble | Branch :: find_na_nd
 
clear; clc; format compact; clf; close all;
 
%% ------------------------------------------------------- input parameters
 
idPHFig = 1;
idVFig  = 2; % Feedrate
idKFig  = 3; % Curvature
idAFig  = 4; % Acceleration
idJFig  = 5; % Jerk (Acceleration of Acceleration)
 
scalarFactor = 100;

% ph Initial and Final points
pI = [1 1]*scalarFactor;   %iniziali sono 1,1_5,2_1,1_1,1
pF = [3 1]*scalarFactor;
 
% ph initial and final tangents
T0 = [1 -1]*scalarFactor; %*30
T1 = [-1 -1]*scalarFactor;  %*20
%T0 = [0.44 -0.29]*scalarFactor; %*30
%T1 = [0 -2.9]*scalarFactor;  %*20



% T0 = [-1 2]*30*scalarFactor;
% T1 = [1 1]*20*scalarFactor;

% [1 1] [1 1]
% [0 1] [1 1]
% [0 1] [1 1]
% [-1 -1] [-1 -1]

% ph solution
useSol_1    = false;
 
% Feedrate values (mm/sec)
Vm = 250;                 % commanded feedrate (m-antain)

VAi = 0;   VAf = Vm;    % Acceleration feedrate init/final
dVAi = 0;  dVAf = 0;    % first derivatives Acceleration feedrate init/final
ddVAi = 0; ddVAf = 0;   % second derivatives Acceleration feedrate init/final

VDi = Vm;   VDf = 0;    % Deceleration feedrate init/final
dVDi = Vm;  dVDf = 0;   % first derivatives Deceleration feedrate init/final
ddVDi = Vm; ddVDf = 0;  % second derivatives Deceleration feedrate init/final

% sampling time (delta t)
deltaT = 0.02; % sec
Ts = deltaT;
 
% Feedrate profile type
% 0 - C0
% 1 - C1
% 2 - C2
VType = 2;

% Tresholds 
Amax  = 800;     % mm/sec^2 Acceleration
Jmax  = 26000;   % mm/sec^2 Jerk
delta = 1;       % mu m (micrometro) chord-error
Ktrhld = 3;

useSigned_K = false;

% colors
PHColor   = 'k';
VaccColor = 'b';
VdecColor = 'g';
 
nTab = 100;
 
%% -------------------------------------------------------------- set input
 
% coefficients of feedrate profile
cVacc = []; % acceleration
cVdec = []; % deceleration
cVcon = []; % -
switch(VType)
    case 0
        Vdegree = 1;
        cVacc = [cVacc VAi VAf];
        cVdec = [cVdec VDi VDf];
    case 1
        Vdegree = 3;
        cVacc = [cVacc VAi VAi VAf VAf];
        cVdec = [cVdec VDi VDi VDf VDf];
    case 2
        Vdegree = 5;
        cVacc = [cVacc VAi VAi VAi VAf VAf VAf];
        cVdec = [cVdec VDi VDi VDi VDf VDf VDf];
end
 
%% ---------------------------------------------------------------- ph init
 
myPH = phCurve(pI, pF, T0, T1, useSol_1);
 
[xi, Cxi] = myPH.computeCurve(nTab);
 
DL  = myPH.computeDisplacementLength();
AL  = myPH.computeArcLength();
CPL = myPH.computeContrPolyLength();
 
%% -------------------------------------------------- K-peaks finding phase

% -- compute K
[par, K, xiKPeaks, KPeaks] = myPH.computeCurvature(useSigned_K,100);

% -- normalize K
[K, KPeaks] = phCurve.normalizeCurvature(K, KPeaks);

% -- plot normalized K
myPH.plotCurvature(par, K, 'k', xiKPeaks, KPeaks, 'r', idKFig);

%% ------------------------------------------------------ plot PH + K-Peaks

% ---- PH
% myPH.plotDispacement(1,'b');
% myPH.plotControlPolygon(1,'b');
phCurve.plotCurve(Cxi, PHColor, 3, idPHFig);

% ---- PH K-Peak on Curvature
for i=1:size(xiKPeaks,2)
    [xi, Cxi] = myPH.computePointOnCurve(xiKPeaks(i));
    figure(idPHFig); hold on;
    plot(Cxi(1),Cxi(2), 'or');
    hold off;
end
    
%% ---- Step I & II: detection of critical points + bounday feedrate values

[xi_CPList, V_CPList] = phCurve.findBoundaryFeedrateVals(delta, Vm, Ts, Amax, Jmax, KPeaks, xiKPeaks, idKFig);

subSegLengthList = myPH.computeSegmentLengths([0 xi_CPList 1]);

% -- foreach ph segment - initialize to 0 all der and der der 
% AccValue{i} = [dViAcc dVfAcc ddViAcc ddVfAcc TAcc]
% DecValue{i} = [dViDec dVfDec ddViDec ddVfDec TDec]
AccValue{1} = [];
DecValue{1} = [];
for i = 1:size(V_CPList,2)+1
    for j=1:5
        AccValue{i}(j) = 0;
        DecValue{i}(j) = 0;
    end
end

% compute coefficient acc/dec for each curve segment.
% cVacCell contains a cell with n array where n is the number of segment of
% ph curve detected in Step I & II. 
% Each array contains the m+1 Berstein coefficients of the acceleration 
% feedrate profile of degree m of the related ph segment.
[cVacCell, cVdeCell] = phCurve.initCoeffVAccDec(VType, Vm, V_CPList, ...
    [], false);

%% ----------------------- Step III: find sampling number for acc/dec Na Nd

B1_VacCell{1} = [];
B1_VdeCell{1} = [];
B2_VacCell{1} = [];
B2_VdeCell{1} = [];
Na_Nd_Nc_Cell{1} = [];
Sa_Sd_Sc_Cell{1} = [];
Ta_Td_Tc_Cell{1} = [];
% for each ph segment
for i = 1:size(cVacCell,2) 
    
    AccCoeffs = cVacCell{i};
    DecCoeffs = cVdeCell{i};
    B1_VacCell{i} = [];
    B1_VdeCell{i} = [];
    B2_VacCell{i} = [];
    B2_VdeCell{i} = [];
    deg = size(AccCoeffs,2)-1;
    
    % compute the second derivative coefficients (for Bound #1)
    for j = 1:size(AccCoeffs,2)-2
        coeffMolt = size(AccCoeffs,2)-2;
        
        cAcc = coeffMolt * (AccCoeffs(j+2) - 2*AccCoeffs(j+1) + AccCoeffs(j));
        B1_VacCell{i} = [B1_VacCell{i} cAcc];
        
        cDec = coeffMolt * (DecCoeffs(j+2) - 2*DecCoeffs(j+1) + DecCoeffs(j));
        B1_VdeCell{i} = [B1_VdeCell{i} cDec];
    end
    
    % compute the third derivative coefficients (for Bound #2)
    for j = 1:size(AccCoeffs,2)-3
        coeffMolt = size(AccCoeffs,2)-3;
        
        cAcc = coeffMolt * (AccCoeffs(j+3) - 3*AccCoeffs(j+2) + 3*AccCoeffs(j+1) - AccCoeffs(j));
        B2_VacCell{i} = [B2_VacCell{i} cAcc];
        
        cDec = coeffMolt * (DecCoeffs(j+3) - 3*DecCoeffs(j+2) + 3*DecCoeffs(j+1) - DecCoeffs(j));
        B2_VdeCell{i} = [B2_VdeCell{i} cDec];
    end
    
    disp(['>> PH Segment #' num2str(i)]);
    disp(['Acc. Coeff: ' num2str(AccCoeffs)]);
    disp(['B1   Coeff: ' num2str(B1_VacCell{i})]);
    disp(['B2   Coeff: ' num2str(B2_VacCell{i})]);
       
    disp(['Dec. Coeff: ' num2str(DecCoeffs)]);
    disp(['B1   Coeff: ' num2str(B1_VdeCell{i})]);
    disp(['B2   Coeff: ' num2str(B2_VdeCell{i})]);
    
    disp(' ');
    
    % -- Bound #1 on 1st derivative (ACC/DEC)
    
    % find ACC roots i-th ph segment (B1) => max B1
    powCoeff_B1_ACC = phCurve.berstainToPowBasis(B1_VacCell{i});
    rootB1_ACC = roots(powCoeff_B1_ACC);
    maxB1_ACC = -10000000000;
    maxr = 0;
    coeffB1_ACC = diff(AccCoeffs); % compute coeff for maxB1
    for j = 1:size(rootB1_ACC,1)
        evalRoot = phCurve.DeCasteljau(coeffB1_ACC, rootB1_ACC(j));
        evalRoot = abs(evalRoot);
        if(isreal(rootB1_ACC(j))&& rootB1_ACC(j)<1 &&rootB1_ACC(j)>0 && evalRoot > maxB1_ACC)
           maxB1_ACC = evalRoot;
           maxr=rootB1_ACC(j);
        end
    end
     if coeffB1_ACC(1) >maxB1_ACC
         maxB1_ACC = coeffB1_ACC(1);
     end
     if coeffB1_ACC(end) >maxB1_ACC
         maxB1_ACC = coeffB1_ACC(end);
     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tauval = linspace(0,1,501);
    % for k=1:501
    %   Vdaval(k)=phCurve.DeCasteljau(coeffB1_ACC,tauval(k));
    % end
    % figure(9+i)
    % plot(tauval,Vdaval,'m',maxr,maxB1_ACC,'ko')    
    
    % find DEC root i-th ph segment (B1) => max B1
    powCoeff_B1_DEC = phCurve.berstainToPowBasis(B1_VdeCell{i});
    rootB1_DEC = roots(powCoeff_B1_DEC);
    maxB1_DEC = -10000000000;
    coeffB1_DEC = diff(DecCoeffs); % compute coeff for maxB1
    for j = 1:size(rootB1_DEC,1)
       evalRoot = phCurve.DeCasteljau(coeffB1_DEC, rootB1_DEC(j));
       evalRoot = abs(evalRoot);
       if(isreal(rootB1_DEC(j))&& rootB1_DEC(j)<1 &&rootB1_DEC(j)>0 && evalRoot > maxB1_DEC)
          maxB1_DEC = evalRoot; 
       end
    end
     if coeffB1_DEC(1) >maxB1_DEC
         maxB1_DEC = coeffB1_DEC(1);
     end
     if coeffB1_DEC(end) >maxB1_DEC
         maxB1_DEC = coeffB1_DEC(end);
     end
    
    
    % -- Bound #2 on 2nd derivative (ACC/DEC)
    maxB2_ACC = [];
    maxB2_DEC = [];
    if(size(B2_VacCell{i},2) >= 2)
        
        % find ACC roots i-th ph segment (B2) => max B2
        powCoeff_B2_ACC = phCurve.berstainToPowBasis(B2_VacCell{i});
        rootB2_ACC = roots(B2_VacCell{i});
        maxB2_ACC = -10000000000;
        coeffB2_ACC = phCurve.commputeDeltaDeltaCoeff(AccCoeffs);
        for j = 1:size(rootB2_ACC,1)
            evalRoot = phCurve.DeCasteljau(coeffB2_ACC, rootB2_ACC(j));
            evalRoot = abs(evalRoot);
            if(isreal(rootB2_ACC(j))&& rootB2_ACC(j)<1 &&rootB2_ACC(j)>0 && evalRoot > maxB2_ACC)
                maxB2_ACC = evalRoot;
            end
        end
        if coeffB2_ACC(1) >maxB2_ACC
            maxB2_ACC = coeffB2_ACC(1);
        end
        if coeffB2_ACC(end) >maxB2_ACC
            maxB2_ACC = coeffB2_ACC(end);
        end
        
        % find DEC roots i-th ph segment (B2) => max B2
        powCoeff_B2_DEC = phCurve.berstainToPowBasis(B2_VdeCell{i});
        rootB2_DEC = roots(powCoeff_B2_DEC);
        maxB2_DEC = -10000000000;
        coeffB2_DEC = phCurve.commputeDeltaDeltaCoeff(DecCoeffs);
        for j = 1:size(rootB2_DEC,1)
            evalRoot = phCurve.DeCasteljau(coeffB2_DEC, rootB2_DEC(j));
            evalRoot = abs(evalRoot);
            if(isreal(rootB2_DEC(j))&& rootB2_DEC(j)<1 &&rootB2_DEC(j)>0 && evalRoot > maxB2_DEC)
                maxB2_DEC = evalRoot;
            end
        end
        if coeffB2_DEC(1) >maxB2_DEC
            maxB2_DEC = coeffB2_DEC(1);
        end
        if coeffB2_DEC(end) >maxB2_DEC
            maxB2_DEC = coeffB2_DEC(end);
        end
        
    end
    
    % -- Find N_a N_d of i-th ph segment
    Na_B1_ACC = (deg/(Amax*deltaT)) * maxB1_ACC; 
    Na_B2_ACC = [];
    if(size(maxB2_ACC,1) ~= 0)
        Na_B2_ACC = sqrt( (deg*(deg-1))/(deltaT^2 * Jmax) * maxB2_ACC);
    end
    Na = max([Na_B1_ACC Na_B2_ACC]);
    
    Na_B1_DEC = (deg/(Amax*deltaT)) * maxB1_DEC; 
    Na_B2_DEC = [];
    if(size(maxB2_DEC,1) ~= 0)
        Na_B2_DEC = sqrt( (deg*(deg-1))/(deltaT^2 * Jmax) * maxB2_DEC);
    end
    Nd = max([Na_B1_DEC Na_B2_DEC]);
    
    % -- Find N_c of i-th ph segment
    % cVacCell{i} = inistial feedrate V_init of i-th segment
    
    Sa = sum(cVacCell{i})/size(cVacCell{i},2) * Na *Ts;
    Sd = sum(cVdeCell{i})/size(cVdeCell{i},2) * Nd *Ts;
    
    Sc = subSegLengthList(i) - Sa - Sd;
    
    Nc = Sc/(Vm*Ts);
    
    Na_Nd_Nc_Cell{i} = [Na Nd Nc];
    
    Sa_Sd_Sc_Cell{i} = [Sa Sd Sc];
    
    Ta = deltaT * Na;
    Td = deltaT * Nd;
    Tc = deltaT * Nc;
    
    Ta_Td_Tc_Cell{i} = [Ta Td Tc];
    
end%_foreach_ph_segments

%% ------------------------------ Step IIIb: changing feedrate coefficients
% 
% for i = 1:size(cVacCell,2) 
%     
%     switch i
%         case 1
%             
%             % Acceleration (init|final) Derivatives
%             dViAcc  = 0;
%             dVfAcc  = 0;
%             ddViAcc = 0;
%             ddVfAcc = 0;
%             
%             % Deceleration (init|final) Derivatives
%             dViDec  = 0;
%             dVfDec  = 400;
%             ddViDec = 0;
%             ddVfDec = 0;
%             
%         case 2
%             
%             % Acceleration (init|final) Derivatives
%             dViAcc  = 400;
%             dVfAcc  = 0;
%             ddViAcc = 1;
%             ddVfAcc = 1;
%             
%             % Deceleration (init|final) Derivatives
%             dViDec  = 0;
%             dVfDec  = -300;
%             ddViDec = 1;
%             ddVfDec = 1;
%             
%         case 3
%             
%             % Acceleration (init|final) Derivatives
%             dViAcc  = -300;
%             dVfAcc  = 0;
%             ddViAcc = 0;
%             ddVfAcc = 0;
%             
%             % Deceleration (init|final) Derivatives
%             dViDec  = 0;
%             dVfDec  = 0;
%             ddViDec = 0;
%             ddVfDec = 0;
%     end
%     
%     TAcc = Ta_Td_Tc_Cell{i}(1);
%     TDec = Ta_Td_Tc_Cell{i}(2);
%     i
%     AcDer = [dViAcc dVfAcc ddViAcc ddVfAcc TAcc]
%     DeDer = [dViDec dVfDec ddViDec ddVfDec TDec]
%     AccValue{i} = [dViAcc dVfAcc ddViAcc ddVfAcc TAcc];
%     DecValue{i} = [dViDec dVfDec ddViDec ddVfDec TDec];
%     disp '';
%     disp '';
% end    
% 
% 
% for i = 1:size(cVacCell,2)
%     disp(['>> PH Segment #' num2str(i)]);
%     disp(['Acc. Coeff: ' num2str(cVacCell{i})]);
%     disp(['Dec. Coeff: ' num2str(cVdeCell{i})]);
%     disp(' ');
% end
% 
% [cVacCell, cVdeCell] = phCurve.initCoeffVAccDec(VType, Vm, V_CPList, ...
%     AccValue, ...
%     DecValue);
% disp('----');
% for i = 1:size(cVacCell,2)
%     disp(['>> PH Segment #' num2str(i)]);
%     disp(['Acc. Coeff: ' num2str(cVacCell{i})]);
%     disp(['Dec. Coeff: ' num2str(cVdeCell{i})]);
%     disp(' ');
% end
%% ------------------ Step IV: find reference parametric values (tau -> xi)

for i_PHSeg = 1:size(cVacCell,2)
    
    Sacc = Sa_Sd_Sc_Cell{i_PHSeg}(1);
    Sdec = Sa_Sd_Sc_Cell{i_PHSeg}(2);
    Scon = Sa_Sd_Sc_Cell{i_PHSeg}(3);
    
    cVacc = cVacCell{i_PHSeg};
    cVdec = cVdeCell{i_PHSeg};
    
    prevSegLen = sum(subSegLengthList(1:i_PHSeg-1));
    
    Tacc = Ta_Td_Tc_Cell{i_PHSeg}(1);
    Tdec = Ta_Td_Tc_Cell{i_PHSeg}(2);
    Tcon = Ta_Td_Tc_Cell{i_PHSeg}(3);
    
    
    % ---- compute F coefficients for ACC/DEC-eleration phase
    cFacc = phCurve.computeFcoeff(cVacc);
    disp(['coefficients of Facc: ' num2str(cFacc)]);
    
    cFdec = phCurve.computeFcoeff(cVdec);
    disp(['coefficients of Fdec: ' num2str(cFdec)]);
    
    % ---- compute Times (Time acceleration, deceleration, constant) (pag. 635)
    AL = subSegLengthList(i_PHSeg);
%     Tacc = Sacc/cFacc(end);
%     Tdec = Sdec/cFdec(end);
%     Tcon = (AL/Vm) - ((Tacc+Tdec)/2);
    T(i_PHSeg) = Tacc + Tdec + Tcon;
    disp('>> Acc | Con | Dec Times | Total Time')
    disp([num2str(Tacc) ' -> ' num2str(Tcon) ' -> ' num2str(Tdec) ' = ' num2str(T(i_PHSeg))]);
    
    % ---- compute parameter values for feedrate profile
    
    prevT = sum(T(1:i_PHSeg-1));
    
    % -- param. values for plotting
    ta = linspace(prevT, prevT+Tacc, 100); 
    td = linspace(prevT, prevT+Tdec, 100); 
    tc = linspace(prevT, prevT+Tcon, 100); 
    
    % -- param. values for computation
    taN = linspace(0, 1, 100);
    tdN = linspace(0, 1, 100);
    tcN = linspace(0, 1, 100);
    
    % -------------------------------------------------------- ACCeleration
    % -- COMPUTE acceleration profile (V, V', V'')
    [Vta] = phCurve.evalFeedRateProfile(cVacc, taN);
    
    [DVta] = phCurve.evalDerFeedrate(cVacc, taN);
    DVta=DVta/Tacc;
    if VType ~= 0
        [DDVta] = phCurve.evalDerDerFeedrate(cVacc, taN);
        DDVta=DDVta/Tacc^2;
    end
    
    % -- PLOT acceleration profile (V, V', V'')
    phCurve.plotFeedRateProfile(ta, Vta,VaccColor,idVFig);
    phCurve.plotDerFeedrate(ta, DVta,VaccColor,idAFig);
    if VType ~= 0
        phCurve.plotDerFeedrate(ta, DDVta,VaccColor,idJFig);
    end
    
    % -------------------------------------------------------- DECeleration
    % -- COMPUTE deceleration profile (V, V', V'')
    [Vtd] = phCurve.evalFeedRateProfile(cVdec, tdN);
    [DVtd] = phCurve.evalDerFeedrate(cVdec, tdN);
    DVtd = DVtd/Tdec;
    if VType ~= 0
        [DDVtd] = phCurve.evalDerDerFeedrate(cVdec, tdN);
         DDVtd=DDVtd/Tdec^2;
    end
    % -- PLOT deceleration profile (V, V', V'')
    phCurve.plotFeedRateProfile(Tacc+Tcon+td, Vtd,VdecColor,idVFig);
    phCurve.plotDerFeedrate(Tacc+Tcon+td, DVtd,VdecColor, idAFig);
    if VType ~= 0
        phCurve.plotDerFeedrate(Tacc+Tcon+td, DDVtd,VdecColor,idJFig);
    end
    
    % --------------- compute reference parameter values for ACCeleration phase
    
    % number of sampling time for acceleration
    XiKacc = [];
    XiKacc(1) = [0];
    N_Xi = floor(Tacc/deltaT);
    incr = deltaT/Tacc;
    for k = 1:N_Xi+1
        
        % normalized sampling time to evaluate integral F of (29.21)
        tkN = (k-1)*incr;
        
        F_tkN = Tacc * phCurve.DeCasteljau(cFacc,tkN) + prevSegLen;
        
        XiKacc(k+1) = myPH.newtonXi(F_tkN, XiKacc(k), 100, exp(-5));
        tkN;
    end
    XiKacc;
    
    % --------------- compute reference parameter values for ACCeleration phase
    
    % number of sampling time for decelaration
    XiKdec = [];
    XiKdec(1) = XiKacc(end);
    N_Xi = floor(Tdec/deltaT);
    for k = 1:N_Xi+1
        
        % normalized sampling time to evaluate integral F of (29.21)
        tkN = (k-1)*deltaT/Tdec;
        % tkN = ((k-1)*deltat - (Tacc+Tcon))/()
        
        F_tkN = Tdec * phCurve.DeCasteljau(cFdec,tkN) + (Sacc + Scon) + prevSegLen;%(AL/3*2); % <-- MOD
        
        XiKdec(k+1) = myPH.newtonXi(F_tkN, XiKdec(k), 100, exp(-5));
        
    end
    XiKdec(1) = XiKdec(2); %update the effective starting value
    
    XiKdec;
    
    % -- Lengths
    % AccLength = phCurve.DeCasteljau(myPH.s, XiKacc(end)) - phCurve.DeCasteljau(myPH.s, XiKacc(1));
    % ConLength = phCurve.DeCasteljau(myPH.s, XiKdec(2)) - phCurve.DeCasteljau(myPH.s, XiKacc(end));
    % DecLength = phCurve.DeCasteljau(myPH.s, 1) - phCurve.DeCasteljau(myPH.s, XiKdec(2));
    % fprintf ( 1, 'ARC Length:\n - Sacc %3.5f | %3.5f \n - Scon %3.5f | %3.5f\n - Sdec %3.5f | %3.5f\n - Stot %3.5f | %3.5f\n', ...
    %     Sacc, AccLength, Scon, ConLength, Sdec, DecLength, AL, AccLength + DecLength + ConLength);
    % AccLength + DecLength + ConLength
    % AL
    
    % ------------------------------------------------- PH + Interp. Ref Points
        
    % ---- interpolation points
    myPH.plotInterpReferencePoints(XiKacc,'m',idPHFig);
    myPH.plotInterpReferencePoints(XiKdec,'g',idPHFig);
    
end

%% ---------------------------------------------------------- plot settings

% -- graphic settings for PH CURVE figure
figure(idPHFig); hold on;
axis equal; 
grid on; grid minor;
set(gcf, 'units','normalized','outerposition',[0 0.55 0.5 0.4]);
set(gcf, 'Color', [1 1 1]);
title('ph curve');
hold off;

% -- graphic settings for ph CURVATURE figure
figure(idKFig); hold on;
grid on; grid minor;
set(gcf, 'units','normalized','outerposition',[0.55 0.55 0.45 0.4]);
set(gcf, 'Color', [1 1 1]);
title('ph curve curvature');
hold off;
 
% -- graphic settings for FEEDRATE (V) profile figure
Vfig = figure(idVFig); hold on;
% axis equal; 
grid on; grid minor;
set(gcf, 'units','normalized','outerposition',[0 0 0.33 0.5]);
set(gcf, 'Color', [1 1 1]);
title('feedrate (V)');
hold off;
 
% -- graphic settings for ACCELERATION (V') profile figure
figure(idAFig); hold on;
% axis equal; 
grid on; grid minor;
set(gcf, 'units','normalized','outerposition',[0.33 0 0.33 0.5]);
set(gcf, 'Color', [1 1 1])
title('acceleration (dV)');
hold off;

% -- graphic settings for JERK (V'') profile figure
figure(idJFig); hold on;
% axis equal; 
grid on; grid minor;
set(gcf, 'units','normalized','outerposition',[0.66 0 0.33 0.5]);
set(gcf, 'Color', [1 1 1]);
title('jerk (ddV)');
hold off;

% movegui(idPHFig,'northwest');
% movegui(idKFig, 'southwest');
% movegui(idVFig, 'northeast');
% movegui(idAFig, 'east');
% movegui(idJFig, 'southeast');

Tfull=0;
for i=1:size(Ta_Td_Tc_Cell,2)
  Tfull= Tfull+ sum(Ta_Td_Tc_Cell{i});
end
disp(['Tfull = ',num2str(Tfull)])
 