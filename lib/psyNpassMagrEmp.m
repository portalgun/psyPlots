function [PA, PC, PAboot, PCboot] = psyNpassMagrEmp(stdX,cmpX,RcmpChs,nPssInd,Magr,stdIind,cmpIind,nPssInd2use,nBootStrp,bPLOT,CIsz)

% function [PA, PC, PAboot, PCboot] = psyNpassMagrEmp(stdX,cmpX,RcmpChs,nPssInd,Magr,stdIind,cmpIind,nPssInd2use,nBootStrp,bPLOT,CIsz)
%
% computes empirical agreement data in an n-tuple-pass experiment
%
% example call: % DATA FROM BMC (PASSES 1-4)
%               S = loadPSYdataSPDall('JND',['BMC'],'NAT',{[1 3 8 10:16 17:26 27 31:32 34:50]},'server');
%               % DATA FROM JDB (PASSES 1-3)
%               S = loadPSYdataSPDall('JND',['JDB'],'NAT',{6:35},'server');
%
%               % PROPORTION CMP CHOSEN VS. PROPORTION AGREEMENT
%               [PA, PC, PAboot, PCboot] = psyNpassMagrEmp(S.stdX,S.cmpX,S.RcmpChs,S.nPssInd,2,S.stdIind,S.cmpIind,[1 2],500,1,90)
%
% stdX:        std speed indices
% cmpX:        cmp speed indices
% RcmpChs:     boolean vector indicating whether cmp stim chosen
% nPssInd:     vector that tags each trial by the pass it came from
% Magr:        M-agree value
% stdIind:     standard stimulus indices
% cmpIind:     comparison stimulus indices
% nPssInd2use: pass indices to use (e.g. if only want to use passes
%              1 and 3, argument is [1 3])
% nBootStrp   : number of bootstraps
% bPLOT:       plot or not
%              1 -> plot
%              0 -> not
% CIsz:        confidence intervals (in percentage, e.g. 90, 95)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PA:          proportion agreement
% PC:          proportion comparison chosen

% outputs:
%          PA: proportion agreement
%          PC: proportion comparison chosen
% cSubplot:       optional subplot count.  This function plots two subplots. cSubplot indicates the index of the first plot if you are building a subplotting routing with more subplots than this function provides.
% tSubplot:       optional max sublot index. Must be defined if using cSubplot. *NOTE, by assigning a value to this parameter, you will be chaning the plot routines to those found in psyFitDecisionVariableCorr.
% bQuiet,         optional paramter to suppress output

if ~exist('nPssInd2use','var') || isempty(nPssInd2use) nPssInd2use =                              [1 2]; end
if ~exist('stdIind','var')     || isempty(stdIind)     stdIind     = [1:numel(stdX)/2 1:numel(stdX)/2]'; end
if ~exist('cmpIind','var')     || isempty(cmpIind)     cmpIind     = [1:numel(cmpX)/2 1:numel(cmpX)/2]'; end
if ~exist('nBootStrp','var')   || isempty(nBootStrp)   nBootStrp   =                                 []; end
if ~exist('bPLOT','var')       || isempty(bPLOT)       bPLOT       =                                  0; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRELIMINARIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INDICES FOR SELECTING TRIALS ACCORDING TO PASS
ind2structSelect = ismember(nPssInd,nPssInd2use);

% CREATE NEW STRUCT THAT ONLY HAS TRIALS FROM THOSE PASSES
nPssInd = nPssInd(ind2structSelect);
stdIind = stdIind(ind2structSelect);
cmpIind = cmpIind(ind2structSelect);
stdX = stdX(ind2structSelect);
cmpX = cmpX(ind2structSelect);
RcmpChs = RcmpChs(ind2structSelect);

% N IS JUST THE LARGEST PASS NUMBER
N = length(unique(nPssInd));

% UNIQUE PASS INDICES
nPssIndUnq = unique(nPssInd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE SURE TRIALS ARE SAME IN EACH PASS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE ARRAYS FOR STORING MOVIE INDICES
stdIindStore = [];
cmpIindStore = [];

for i = 1:length(nPssIndUnq) % FOR EACH UNIQUE PASS
   % GRAB STD AND CMP MOVIES
   bIndTmp = nPssInd == nPssIndUnq(i); 
   stdIindStore(:,i) = stdIind(bIndTmp);
   cmpIindStore(:,i) = cmpIind(bIndTmp);
end

% INITIALIZE ARRAYS FOR STORING SORTED MOVIE INDICES
stdIindSorted = [];
cmpIindSorted = [];

for i = 1:size(stdIindStore,2) % FOR EACH PASS
   % GRAB STD AND CMP MOVIE INDICES
   stdIindTmp = stdIindStore(:,i);
   cmpIindTmp = cmpIindStore(:,i);
   % SORT STIMULUS PAIRS FOR EACH TRIAL
   [~,alignInd] = sortrows([stdIindTmp cmpIindTmp]);
   % SORT THEM AND STORE
   stdIindSorted(:,i) = stdIindTmp(alignInd);
   cmpIindSorted(:,i) = cmpIindTmp(alignInd);
end

% CHECK THAT SAME STD AND CMP MOVIES WERE USED IN EACH PASS
if ~isequal(sum(stdIindSorted,2)./size(stdIindSorted,2),stdIindSorted(:,1)) || ...
   ~isequal(sum(cmpIindSorted,2)./size(cmpIindSorted,2),cmpIindSorted(:,1))
   error('psyNpassMagrEmp: Stimuli not same across passes!');
end

%%%%%%%%%%%%%%%%%%%%
% COMPUTE AGREEMENT
%%%%%%%%%%%%%%%%%%%%

% Note from Ben: these loops might not be the easiest to read (I have tried), but the
% logic is straightforward. First, I loop over all of the standard speeds
% and grab the relevant fields in the struct for each loop. This is the 'i'
% loop, and all the relevant fields are labeled 'tmp' for temporary. Within
% each of those standard speed ('i') loops, I loop over all comparison speeds and
% grab the relevant fields. This is the 'j' loop, and all of the relevant
% fields are now labeled 'tmpTWO' to distinguish them from the 'tmp' fields in the
% previous (higher level) loop. Then, in the 'k' loop, I loop over passes
% and fields within thise loop are labeled 'tmpTHREE' to distinguish them
% from fields in the previous 'j' loop.

% UNIQUE STD SPEEDS
stdXunq = unique(abs(stdX)); 

PA = [];
PC = [];
PAboot = [];
PCboot = [];

for i = 1:length(stdXunq) % FOR EACH UNIQUE STD SPEED
   % GRAB INDICES WITH THAT STD SPEED
   bIndTmp = abs(abs(stdX)-stdXunq(i)) < 0.01;
   % GRAB CMP X VECTOR
   cmpXtmp =    cmpX(bIndTmp);
   % GRAB CMP CHS VECTOR
   cmpChsTmp =  RcmpChs(bIndTmp);
   % GRAB PASS INDEX VECTOR
   nPssIndTmp = nPssInd(bIndTmp);
   % GRAB MOVIE INDICES
   stdIindTmp = stdIind(bIndTmp);
   cmpIindTmp = cmpIind(bIndTmp);
   % UNIQUE CMP CHS
   cmpXunqTmp = unique(abs(cmpXtmp));
   for j = 1:length(cmpXunqTmp) % FOR EACH UNIQUE CMP SPEED
      % GRAB INDICES WITH THAT SPEED
      bIndTmpTWO = abs(abs(cmpXtmp)-cmpXunqTmp(j)) < 0.01;
      % GRAB CMP CHS VECTOR
      cmpChsTmpTWO  = cmpChsTmp (bIndTmpTWO);
      % GRAB PASS INDEX VECTOR
      nPssIndTmpTWO = nPssIndTmp(bIndTmpTWO);
      % GRAB MOVIE INDEXES
      stdIindTmpTWO = stdIindTmp(bIndTmpTWO);
      cmpIindTmpTWO = cmpIindTmp(bIndTmpTWO);
      % UNIQUE PASS INDICES
      nPssIndTmpUnqTWO = unique(nPssIndTmpTWO);
      % INITIALIZE ARRAY FOR STORING CMP CHS VALUES. THIS WILL BE A
      % numTrials X N MATRIX--EACH COLUMN CONTAINS THE CMP CHS VALUES FOR
      % A SINGLE PASS
      cmpChsPassKtmpTWO = [];
      for k = 1:length(nPssIndTmpUnqTWO) % FOR EACH UNIQUE MOVIE PASS
         % SPLIT STD AND CMP MOVIE INDICES ACCORDING TO THEIR PASSES
         stdIindPassKtmpTHREE = stdIindTmpTWO(nPssIndTmpTWO==nPssIndTmpUnqTWO(k));
         cmpIindPassKtmpTHREE = cmpIindTmpTWO(nPssIndTmpTWO==nPssIndTmpUnqTWO(k));
         % ALIGN ACCORDING TO STD AND CMP SPEEDS
         [~,alignInd] = sortrows([stdIindPassKtmpTHREE cmpIindPassKtmpTHREE]);
         % SPLIT INTO MULTIPLE PASSES
         cmpChsPassKtmpTHREE = cmpChsTmpTWO(nPssIndTmpTWO==nPssIndTmpUnqTWO(k));
         % ALIGN CMP CHS VECTORS BY TRIAL IDENTITY
         cmpChsPassKtmpTHREE = cmpChsPassKtmpTHREE(alignInd);
         % STORE IN MATRIX
         cmpChsPassKtmpTWO(:,k) = cmpChsPassKtmpTHREE;
      end
      % COMPUTE AGREEMENT AND PROP CMP CHS
      PA(i,j) = (sum(sum(cmpChsPassKtmpTWO==0,2)==Magr) + sum(sum(cmpChsPassKtmpTWO==1,2)==Magr))./length(alignInd);
      PC(i,j) = sum(cmpChsTmpTWO)./length(cmpChsTmpTWO);
      if ~isempty(nBootStrp)
          for m = 1:nBootStrp % FOR NUMBER OF BOOTSTRAPS
             % BOOTSTRAP P(C)
             cmpChsTmpTWOboot = randsample(cmpChsTmpTWO,length(cmpChsTmpTWO),1);
             PCboot(i,j,m) = sum(cmpChsTmpTWOboot)./length(cmpChsTmpTWOboot);
             % BOOTSTRAP P(A)
             cmpChsPassKtmpTWOboot = cmpChsPassKtmpTWO(randsample(1:size(cmpChsPassKtmpTWO,1),size(cmpChsPassKtmpTWO,1),1),:);
             PAboot(i,j,m) = (sum(sum(cmpChsPassKtmpTWOboot==0,2)==Magr) + sum(sum(cmpChsPassKtmpTWOboot==1,2)==Magr))./length(alignInd);
          end
      end
   end
end

if bPLOT == 1
    % VARIABLE FOR CONVENIENT SCALING AXES (E.G. PROPORTION OR %)
    yScale = 1;
    % NUMBER OF BINOMIAL SIMULATIONS
    numSim = 20000;
    % QUANTILES FOR BOOTSTRAPPED CONFIDENCE INTERVALS
    CI = [0.5.*(1-CIsz/100) 1-0.5.*(1-CIsz/100)];
    % GENERATE BINOMIAL PREDICTION
    [PAbin, PCbin, PAlb, PAub] = psyNpassBinomPred(N,Magr,repmat([-2:0.2:2]',[1 N]),[1 0; 0 1],CI,numSim,length(cmpChsPassKtmpTHREE),0);
    % RESOLUTION OF INTERPOLATION
    PCbinRes = 1./numel(cmpChsPassKtmpTWO);
    % PC VALUES FOR INTERPOLATION
    PCbinInterp = [PCbinRes:PCbinRes:1-PCbinRes]';
    % INTERPOLATE AGREEMENT
    PAbinInterp = spline(PCbin,PAbin,PCbinInterp);
    % INTERPOLATE LOWER BOUND ON AGREEMENT
    PAlbInterp = spline(PCbin,PAlb,PCbinInterp);
    % INTERPOLATE UPPER BOUND ON AGREEMENT
    PAubInterp = spline(PCbin,PAub,PCbinInterp);

    %%%%%%%%%%%%%%
    % PLOT STUFF %
    %%%%%%%%%%%%%%
    figure;
    set(gcf,'Position',[569 518 1501 689]);
    for i = 1:length(stdXunq)
       subplot(2,3,i)
       fill([yScale*PCbinInterp; yScale*flipud(PCbinInterp)],[yScale*PAlbInterp; yScale*PAubInterp],0.9.*[1 1 1],'EdgeColor','none'); hold on;
       plot(yScale*PC(i,:),yScale*PA(i,:),'ks','LineWidth',2,'MarkerSize',15); hold on;
       plot(yScale*PCbinInterp,yScale*PAbinInterp,'--k','LineWidth',1.5);
       axis square;
       Fig.format('% Comparison Chosen','% Agreement',['stdX = ' num2str(round(stdXunq(i),3,'significant'))]);
       set(gca,'LineWidth',1.5);
       ylim([0.4*yScale yScale]);
       set(gca,'XTick',[0:0.25*yScale:yScale]);
    end
elseif bPLOT == 1 && bPlotManual == 1 %IF THE SUBPLOTS FROM THIS FUNCTION ARE COMBINED WITH SUBPLOTS FROM OTHER PLOTTING ROUTINES
    subplot(ceil(tSubplot/2),2,cSubplot)
    fill([yScale*PCbinInterp; yScale*flipud(PCbinInterp)],[yScale*PAlbInterp; yScale*PAubInterp],0.9.*[1 1 1],'EdgeColor','none'); hold on;
    plot(yScale*PC(i,:),yScale*PA(i,:),'ks','LineWidth',2,'MarkerSize',15); hold on;
    plot(yScale*PCbinInterp,yScale*PAbinInterp,'--k','LineWidth',1.5);
    axis square;
    Fig.format('% Comparison Chosen','% Agreement',['stdX = ' num2str(round(stdXunq(i),3,'significant'))]);
    set(gca,'LineWidth',1.5);
    ylim([0.4*yScale yScale]);
    set(gca,'XTick',[0:0.25*yScale:yScale]);
end

% TRANSPOSE
PA = PA';
PC = PC';
PAboot = permute(PAboot,[2 1 3]);
PCboot = permute(PCboot,[2 1 3]);

end

