classdef phCurve
    
    properties
        
        U       = [];
        V       = [];
        W       = [];
        sigma   = [];
        s       = [];
        L       = [];
        CPcmplx = [];
        CP      = [];   % 2 x #ControlPoint
        
        Kdata = [];
        
    end
    
    methods
        
        % -- Constructor
        
        function this = phCurve(pI, pF, T0, T1, useSol_1)
            
            % ---- compute Ph Quintic by hermite interpolation
            
            d0 = complex(T0(1), T0(2));
            d1 = complex(T1(1), T1(2));
            
            pI = complex(pI(1), pI(2));
            pF = complex(pF(1), pF(2));
            
            % -- pag. 530 (25.18)
            w0 = sqrt(d0);
            w2 = sqrt(d1);
            
            coeff_1 = 2/3;
            coeff_2 = w0 + w2;
            coeff_3 = w0^2 + w2^2 + (w2*w0)/3 - 5*(pF-pI);
            
            w1 = roots([coeff_1 coeff_2 coeff_3]);
            
            if(useSol_1)
                w1 = w1(1);
            else
                w1 = w1(2);
            end
            
            % -- pag. 415 (19.12)
            p1 = pI + 1/5 * w0^2;
            p2 = p1 + 1/5 * w0*w1;
            p3 = p2 + 1/15 * (2*w1^2 + w0*w2);
            p4 = p3 + 1/5 * w1*w2;
            p5 = p4 + 1/5 * w2^2;
            CP = [pI p1 p2 p3 p4 p5];
            this.CPcmplx = CP;
            
            this.CP = [real(this.CPcmplx);
                imag(this.CPcmplx)];
            
            % -- save w,u,v coefficents (building polynomials W complex, and coeffs not complex)
            W = [w0 w1 w2];
            U = [real(w0) real(w1) real(w2)];
            V = [imag(w0) imag(w1) imag(w2)];
            this.W = W;
            this.U = U;
            this.V = V;
            
            % -- save sigmas
            this.sigma(1) = U(1)^2 + V(1)^2;
            this.sigma(2) = U(1)*U(2) + V(1)*V(2);
            this.sigma(3) = 2/3* (U(2)^2 + V(2)^2) + 1/3*(U(1)*U(3) + V(1)*V(3));
            this.sigma(4) = U(2)*U(3) + V(2)*V(3);
            this.sigma(5) = U(3)^2 + V(3)^2;
            
            % -- save arc length coefficients s
            this.s(1) = 0;
            this.s(2) = 1/5 * (this.sigma(1));
            this.s(3) = 1/5 * (this.sigma(1) + this.sigma(2));
            this.s(4) = 1/5 * (this.sigma(1) + this.sigma(2) + this.sigma(3));
            this.s(5) = 1/5 * (this.sigma(1) + this.sigma(2) + this.sigma(3) + this.sigma(4));
            this.s(6) = 1/5 * (this.sigma(1) + this.sigma(2) + this.sigma(3) + this.sigma(4) + this.sigma(5));
            
        end
        
        % -- computation function
        
        function xiNext = newtonXi(this, F_kt, xiCurr, itMax, tolx)
            
            fx  = phCurve.DeCasteljau(this.s,     xiCurr) - F_kt;
            f1x = phCurve.DeCasteljau(this.sigma, xiCurr);
            xiNext = xiCurr - (fx / f1x);
            
            i = 0;
            
            while(i < itMax && abs(xiNext - xiCurr) > tolx)
                
                i = i+1;
                xiCurr = xiNext;
                
                fx  = phCurve.DeCasteljau(this.s,     xiCurr) - F_kt;
                f1x = phCurve.DeCasteljau(this.sigma, xiCurr);
                xiNext = xiCurr - (fx / f1x);
                
            end
            
            if abs(xiNext - xiCurr) > tolx
                error(['>> no convergence -> xiCurr: ' num2str(xiCurr) ' | xiNext: ' num2str(xiNext) '\n']);
            end
            
        end
        
        function [xi, Cxi] = computeCurve(this, nTab)
            
            xi = linspace(0,1, nTab);
            Cxi = [];
            
            for i = 1:nTab
                Cxi(:,i) = phCurve.DeCasteljau(this.CP,xi(i));
            end
            
        end
        
        function [xi, Cxi] = computePointOnCurve(this, xi)
            Cxi = phCurve.DeCasteljau(this.CP,xi);
        end
        
        function CPL = computeContrPolyLength(this)
            CPL = 0;
            for i=1:size(this.CP,2)-1
                CPL = CPL + norm(this.CP(:,i+1) - this.CP(:,i));
            end
        end
        
        function DL = computeDisplacementLength(this)
            DL = norm(this.CP(:,end) - this.CP(:,1));
        end
        
        function AL = computeArcLength(this)
            AL = sum(this.sigma)/(size(this.CP,2)-1);
        end
        
        function [xi, K, xiMax, Kmax] = computeCurvature(this, isSigned,nTab)
            
            xi = linspace(0,1,nTab);
            xiMax = [];
            K = [];
            Kmax = [];
            
            % -- implementation of pag.388 (17.16)
            deltaU = 2*[this.U(2)-this.U(1) this.U(3)-this.U(2)];
            deltaV = 2*[this.V(2)-this.V(1) this.V(3)-this.V(2)];
            for i = 1:nTab
                
                Ut = phCurve.DeCasteljau(this.U, xi(i));
                Vt = phCurve.DeCasteljau(this.V, xi(i));
                U1t = phCurve.DeCasteljau(deltaU, xi(i));
                V1t = phCurve.DeCasteljau(deltaV, xi(i));
                sigma2 = phCurve.DeCasteljau(this.sigma, xi(i))^2;
                
                K(i) = 2 * (Ut*V1t - U1t*Vt)/(sigma2);
                
                if(~isSigned)
                    K(i) = abs(K(i));
                    % -- curvature peaks detection
                    if( i>3 && K(i-2)<K(i-1) && K(i-1)>K(i))
                        xiMax = [xiMax xi(i-1)];
                        Kmax = [Kmax K(i-1)];
                    end
                end
            end
            
        end
        
        function [xi, T] = computeParametricSpeed(this, nTab)
            
        end
        
        % -- graphic function
        
        function plotControlPolygon(this, idFig, colorCode)
            figure(idFig); hold on;
            plot(this.CP(1,:), this.CP(2,:),['-o' colorCode]);
            hold off;
        end
        
        function plotDispacement(this, idFig, colorCod)
            figure(idFig); hold on;
            plot([this.CP(1,1) this.CP(1,end)], [this.CP(2,1) this.CP(2,end)], ['-o' colorCod]);
            hold off;
            
        end
        
        function plotInterpReferencePoints(this, Xi, colCod, idFig)
            
            figure(idFig); hold on;
            for i=1:size(Xi,2)
                
                Cxi = phCurve.DeCasteljau(this.CP, Xi(i));
                plot(Cxi(1), Cxi(2),['.' colCod],'MarkerSize',8);
                
            end
            hold off;
            
        end
        
        function plotCurvature(this, xi, K, KColCode, xiMax, KMax, peaksColCode, idFig)
            
            figure(idFig); hold on;
            plot(xi, K, ['-' KColCode]);
            plot(xiMax, KMax, ['o' peaksColCode]);
            hold off;
            
        end
        
        function segArcLen = computeSegmentArcLength(this, xiStart, xiEnd)
            LStart  = phCurve.DeCasteljau(this.s, xiStart);
            LEnd    = phCurve.DeCasteljau(this.s, xiEnd);
            segArcLen = LEnd - LStart;
        end
        
        function subSegLengthList =  computeSegmentLengths(this, xi_CPList)
            subSegLengthList = [];
            for i = 1:size(xi_CPList,2)-1
                subSegLengthList(i) = this.computeSegmentArcLength(xi_CPList(i), xi_CPList(i+1));
            end
        end
        
        
    end
    
    methods (Static)
        
        %
        function Cu = DeCasteljau(P,u)
            
            n = size(P,2);
            Q = [];
            for i= 1:n
                Q(:,i) = P(:,i);
            end
            
            for k = 1:n
                for i = 1: n-k
                    Q(:,i) = (1-u)*Q(:,i) + u*Q(:,i+1);
                end
            end
            
            Cu = Q(:,1);
        end
        
        function plotCurve(Cxi, colorCode, thick,idFig)
            figure(idFig); hold on;
            plot(Cxi(1,:), Cxi(2,:),['-' colorCode], 'LineWidth', thick);
            hold off;
        end
        
        function [evalVti] = evalFeedRateProfile(V,ti)
            
            evalVti = zeros(1,size(ti,2));
            for i=1:size(ti,2)
                Vti = phCurve.DeCasteljau(V,ti(i));
                evalVti(i) = Vti;
            end
            
        end
        
        function [evalV1ti] = evalDerFeedrate(V,ti)
            
            n = size(V,2)-1;
            deltaV = [];
            for i = 1:n
                deltaV(i) = V(i+1) - V(i);
            end
            deltaV = n*deltaV;
            
            evalV1ti = [];
            for i = 1:size(ti,2)
                evalV1ti(i) = phCurve.DeCasteljau(deltaV, ti(i));
            end
            
        end
        
        % compute delta(coefficient vector) = diff(V)
        % V      - [1xn]
        % DeltaV - [1x(n-1)]
        function deltaV = computeDeltaCoeff(V)
            n = size(V,2)-1;
            deltaV = [];
            for i = 1:n
                deltaV(i) = V(i+1) - V(i);
            end
        end
        
        % compute DeltaDelta(coefficient vector)
        % V      - [1xn]
        % DeltaV - [1x(n-2)]
        function deldelV = commputeDeltaDeltaCoeff(V)
            n = size(V,2)-2;
            deldelV = [];
            for i = 1:n
                deldelV(i) = V(i+2) - 2*V(i+1) + V(i);
            end
        end
        
        function [evalV2ti] = evalDerDerFeedrate(V,ti)
            
            n = size(V,2)-2;
            deldelV = [];
            for i = 1:n
                deldelV(i) = V(i+2) - 2*V(i+1) + V(i);
            end
            deldelV = n*(n-1)*deldelV;
            
            evalV2ti = [];
            for i = 1:size(ti,2)
                evalV2ti(i) = phCurve.DeCasteljau(deldelV, ti(i));
            end
            
        end
        
        function plotFeedRateProfile(ti,Vti, colCode, idFig)
            figure(idFig); hold on;
            plot(ti, Vti, ['-' colCode]);
            hold off;
        end
        
        function plotDerFeedrate(ti, DVti, colCode, idFig)
            figure(idFig); hold on;
            plot(ti, DVti, ['-' colCode]);
            hold off;
        end
        
        function [coeffF] = computeFcoeff(coeffV)
            coeffF = [];
            coeffF(1) = 0;
            FDeg = size(coeffV,2);
            FNumCoeff = size(coeffV,2)+1;
            for i = 2:FNumCoeff
                coeffF(i) = 0;
                for j = 1:(i-1)
                    coeffF(i) = coeffF(i) + coeffV(j);
                end
                coeffF(i) = coeffF(i)/FDeg;
            end
        end
        
        function [K, KPeaks] = normalizeCurvature(K, KPeaks)
            K = K/ceil(max(KPeaks));
            KPeaks = KPeaks/ceil(max(KPeaks));
        end
        
        function [ Na_Nd_Nc_Cell, Ta_Td_Tc_Cell, Sa_Sd_Sc_Cell, cVacCell, cVdeCell ] = ...
                findAccDecTimes(cVacCell, cVdeCell, Vm, Amax, Jmax, deltaT, Ts, subSegLengthList)
            
            B1_VacCell{1} = [];
            B1_VdeCell{1} = [];
            B2_VacCell{1} = [];
            B2_VdeCell{1} = [];
            Na_Nd_Nc_Cell{1} = [];
            Sa_Sd_Sc_Cell{1} = [];
            Ta_Td_Tc_Cell{1} = [];
            
            NumVmRefinements = 100;
            VmValueList = linspace(Vm, 0, NumVmRefinements);
            
            
            % for each ACC/DEC Phase
            for i = 1:size(cVacCell,2)
                
                % If there is need to recompute all times due to Nc<0
                for k = 1:NumVmRefinements-1
                    
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
                    
                    if(0)
                        disp(['>> PH Segment #' num2str(i)]);
                        disp(['Acc. Coeff: ' num2str(AccCoeffs)]);
                        disp(['B1   Coeff: ' num2str(B1_VacCell{i})]);
                        disp(['B2   Coeff: ' num2str(B2_VacCell{i})]);
                        
                        disp(['Dec. Coeff: ' num2str(DecCoeffs)]);
                        disp(['B1   Coeff: ' num2str(B1_VdeCell{i})]);
                        disp(['B2   Coeff: ' num2str(B2_VdeCell{i})]);
                        
                        disp(' ');
                    end
                    
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
                    
                    % Manage case there is no constant phase time
                    % if(0)
                    if(Nc < 0)
                        
                        % --- (START BLOCK IGNORED)
                        if(0)% find_intersection_wrong_method
                            
                            % find intersection between feedrate profile
                            % Newton applied on f(t) = Vacc(t) - Vdec(t)
                            % to find roots of f(t) = 0
                            
                            t0 = Ta+Tc;
                            tCurr = t0;
                            
                            % -- plot ACC
                            taN = linspace(0, 1, 100);
                            cVAcc = cVacCell{i};
                            [Vta] = phCurve.evalFeedRateProfile(cVAcc, taN);
                            ta = linspace(0,Ta,100);
                            if(i==1)
                                % phCurve.plotFeedRateProfile(ta, Vta,'g',100);
                                phCurve.plotFeedRateProfile(taN, Vta,'g',100);
                                figure(100); grid on; grid minor;
                            else
                                phCurve.plotFeedRateProfile(taN, Vta,'r',101);
                                figure(101); grid on; grid minor;
                            end
                            
                            cVDec = cVdeCell{i};
                            
                            % -- Newton
                            tol = exp(-9);
                            itMax = 100;
                            
                            fxAcc  = phCurve.DeCasteljau(cVAcc, tCurr/Ta);
                            fxDec  = phCurve.DeCasteljau(cVDec, (tCurr - t0)/(Td));
                            
                            diffCVacc = (cVAcc(2:end) - cVAcc(1:end-1))/Ta;
                            diffCVdec = (cVDec(2:end) - cVDec(1:end-1))/(Td);
                            
                            f1xAcc = 5*phCurve.DeCasteljau(diffCVacc, tCurr/Ta);
                            f1xDec = 5*phCurve.DeCasteljau(diffCVdec, (tCurr - t0)/(Td));
                            
                            tNext = tCurr - ( (fxAcc - fxDec)/(f1xAcc - f1xDec) );
                            
                            j = 0;
                            while(j < itMax && abs(tNext - tCurr) > tol)
                                j = j+1;
                                tCurr = tNext;
                                
                                fxAcc  = phCurve.DeCasteljau(cVAcc, tCurr/Ta);
                                fxDec  = phCurve.DeCasteljau(cVDec, (tCurr - t0)/(Td));
                                
                                diffCVacc = (cVAcc(2:end) - cVAcc(1:end-1))/Ta;
                                diffCVdec = (cVDec(2:end) - cVDec(1:end-1))/(Td);
                                
                                f1xAcc = 5*phCurve.DeCasteljau(diffCVacc, tCurr/Ta);
                                f1xDec = 5*phCurve.DeCasteljau(diffCVdec, (tCurr - t0)/(Td));
                                
                                tNext = tCurr - ( (fxAcc - fxDec)/(f1xAcc - f1xDec) );
                            end
                            
                            if(abs(tNext - tCurr) > tol)
                                warning(['>> no convergence -> tauCurr: ' num2str(tCurr) ' | tauNext: ' num2str(tNext) '\n']);
                            end
                            % -- end Newton
                            
                            tDis = tNext;
                            if(i>1)
                                for k = 1:i-1
                                    tDis =tDis + sum(Ta_Td_Tc_Cell{k});
                                end
                            end
                            figure(3); hold on; plot(tDis, phCurve.DeCasteljau(cVAcc,tNext/Ta),'or');
                            %
                            if(1) % Modifico coefficienti feedrate e ricalcolo i tempi
                                
                                % Vfs = phCurve.DeCasteljau(cVAcc,tNext/Ta);
                                Vfs = 190;
                                Vi = AccCoeffs(1);
                                Vf = DecCoeffs(end);
                                % -- change the constant feedrate value Vm and
                                % compute again Acc/Dec Times
                                % (in realt conoscendo TNext gi conosco Ta Td)
                                cVa3 = alphaAcc(i)*Vi + ((1 - alphaAcc(i)) * Vfs);
                                cDe4 = alphaDec(i)*Vfs + ((1 - alphaDec(i)) * Vf);
                                
                                cVacCell{i} = [Vi  Vi  cVa3 Vfs Vfs Vfs];
                                cVdeCell{i} = [Vfs Vfs Vfs  cDe4 Vf Vf];
                                
                            else % Non modifico i coefficient Feedrate
                                break;
                            end
                            
                        end %_end_find_intersection_wrong_method
                        % --- (END BLOCK IGNORED)
                          
                        % -- further Vm refinements are required
                        
                        %break;
                        
                        % Vfs = phCurve.DeCasteljau(cVAcc,tNext/Ta);
                        Vfs = VmValueList(k+1);
                        %Vi = AccCoeffs(1);
                        %Vf = DecCoeffs(end);
                        
                        % -- change the constant feedrate value Vm and
                        % compute again Acc/Dec Times
                        % (in realt conoscendo TNext gi conosco Ta Td)
                        %cVa3 = alphaAcc(i)*Vi + ((1 - alphaAcc(i)) * Vfs);
                        %cDe4 = alphaDec(i)*Vfs + ((1 - alphaDec(i)) * Vf);
                        
                        cVacCell{i}(4:6) = [Vfs Vfs Vfs];
                        cVdeCell{i}(1:3) = [Vfs Vfs Vfs];
                        
                    else
                        
                        % -- no further Vm refinements are required
                        break;
                        
                    end%_end_if_Nc<0
                    
                end%_for_1:NumVmRefinements
                
            end%_foreach_AccDec_phase
            
            
            
        end
        
        
        
        % Compute the feedrate values on critical points
        %
        % Input:
        % delta     - [double] chord error (default 1)
        % Vm        - [double] commanded feedrate
        % Ts        - [double] sampling time
        % Amax      - [double] maximum acceleration (mm/sec^2)
        % Jmax      - [double] maximum jerk (mm/sec^2)
        % KPeaks    - [1 x #Peaks] array of peaks of curvature
        % xiKPeaks  - [1 x #Peaks] array of parameter value xi -> KPeaks(i)
        % idKFig    - [int] curvature figure id value
        %
        % Output:
        % V_CPList  - [1 x #CriticalPoints] array of feedrate values at
        %             detected critical points
        function [xi_CPList, V_CPList] = findBoundaryFeedrateVals(delta, Vm, Ts, Amax, Jmax, KPeaks, xiKPeaks, idKFig)
            % soglie
            K_thr1 = 8*delta /((Vm*Ts)^2 + 4*delta^2);
            K_thr2 = Amax/Vm^2;
            K_thr3 = sqrt(Jmax/Vm^3);
            Ktrhld = min([K_thr1, K_thr2, K_thr3]);
            figure(idKFig); hold on;
            plot([0,1],[Ktrhld Ktrhld],'-.r');
            hold off;
            
            xi_CPList = []; % list of critical point xi values
            K_CPList  = []; % list of critical point K-urvature values
            V_CPList  = []; % list of critical point feedrate values
            for i = 1:size(KPeaks,2)
                if(KPeaks(i) > Ktrhld)
                    xi_CPList = [xi_CPList xiKPeaks(i)];
                    K_CPList = [K_CPList KPeaks(i)];
                    
                    % find the current critical point feedrate value (frontier feedrate)
                    V_CP_1 = 2/Ts * sqrt( (1/KPeaks(i)^2) - (1/KPeaks(i) - delta)^2 );
                    V_CP_2 = sqrt(Amax/KPeaks(i));
                    V_CP_3 = (Jmax/KPeaks(i)^2)^(1/3);
                    V_CP = min([V_CP_1  V_CP_2  V_CP_3]);
                    V_CPList = [V_CPList V_CP];
                end
            end
        end
        
        
        %
        function [cVacCell, cVdeCell, alphaAcc, alphaDec] = initCoeffVAccDec(VType, Vm, V_CPList, alpha, useAlphaOpt)
            
            if(VType == 2)
                n = length(alpha);
                d = n/2;
                alphaAcc(2:d+1) = alpha(1 : d);
                alphaDec = alpha(d+1 : n);
                alphaDec(d+1) = 0;
            end
            
            % coefficients of feedrate profile
            cVacCell{1} = []; % acceleration
            cVdeCell{1} = []; % deceleration
            
            % for each curve segment
            nSeg = size(V_CPList,2)+1;
            for i = 1:nSeg
                
                switch i
                    case 1 % primo tratto
                        VAi = 0;
                        VAf = Vm;
                        VDi = Vm;
                        VDf = V_CPList(1);
                        
                    case nSeg % ultimo tratto
                        VAi = V_CPList(end);
                        VAf = Vm;
                        VDi = Vm;
                        VDf = 0;
                        
                    otherwise % altri tratti
                        VAi = V_CPList(i-1);
                        VAf = Vm;
                        VDi = Vm;
                        VDf = V_CPList(i);
                end
                
                switch(VType)
                    
                    case 0 % feedrate grado 1
                        cVacCell{i} = [VAi VAf];
                        cVdeCell{i} = [VDi VDf];
                        
                    case 1 % feedrate grado 3
                        cVacCell{i} = [VAi VAi VAf VAf];
                        
                        cVdeCell{i} = [VDi VDi VDf VDf];
                        
                    case 2 % feedrate grado 5
                        
                        if(useAlphaOpt)
                            
                            cVa3 = alphaAcc(i)*VAi + ((1 - alphaAcc(i)) * VAf);
                            cDe4 = alphaDec(i)*VDi + ((1 - alphaDec(i)) * VDf);
                            
                            % without optimization process => JERK DISCONTINUE
                            cVacCell{i} = [VAi VAi cVa3  VAf VAf VAf];
                            cVdeCell{i} = [VDi VDi VDi   cDe4 VDf VDf];
                            
                        else
                            
                            cVacCell{i} = [VAi VAi VAi  VAf VAf VAf];
                            cVdeCell{i} = [VDi VDi VDi  VDf VDf VDf];
                        end
                        
                end
            end
            
        end%_initCoeffVAccDec
        
        % Implementation of the last formulae in http://cagd.cs.byu.edu/~557/text/ch3.pdf
        % Sedemberg CAGD Course note 2016. Convert the arry of berstein
        % coefficient to power base coefficents.
        %
        %Input:
        %   barCoeff    - [1x#Coeff] array of coefficients in Berstain base
        %
        %Output:
        %   powCoeff    - [1x#Coeff] arry of power base coefficients (greate -> lower)
        %
        function [powCoeff] = berstainToPowBasis(berCoeff)
            % checked with: http://cagd.cs.byu.edu/~557/text/ch3.pdf
            powCoeff = [];
            deg = size(berCoeff,2)-1;
            
            for i=1:size(berCoeff,2)
                currCoeff = 0;
                for k = 1:i
                    currCoeff = currCoeff + berCoeff(k)*nchoosek(deg,i-1)*nchoosek(i-1,k-1)*(-1)^(i-k);
                end
                powCoeff(i) = currCoeff;
            end
            
            powCoeff = flip(powCoeff);
        end
        
        function plotConFeedrate(tacc, tdec, Vtacc, Vtdec, dVtacc, dVtdec, ddVtacc, ddVtdec,...
                idFigV, idFigA, idFigJ, colorCode)
            
            figure(idFigV);
            hold on;
            plot([tacc tdec], [Vtacc Vtdec], ['-' colorCode]);
            hold off;
            
            figure(idFigA);
            hold on;
            plot([tacc tdec], [dVtacc dVtdec], ['-' colorCode]);
            hold off;
            
            if(size(ddVtacc,1) ~= 0)
                figure(idFigJ);
                hold on;
                plot([tacc tdec], [ddVtacc ddVtdec], ['-' colorCode]);
                hold off;
            end
            
        end
        
        
    end%_Static
end