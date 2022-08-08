% **************************************************************************************
% This code implements adaptive-bandwidth kernel density estimation and is intended for 
% use on high-throughput sequencing data. It accompanies the paper "Adaptive Bandwidth
% Kernel Density Estimation for Next-Generation Sequencing Data" by P. Ramachandran and
% T. J. Perkins (submitted, 2013). This code is free to use for non-commercial purposes.
%
% USAGE:
% ---------
% Basically, this software takes in a BED file containing the read locations of an
% NGS dataset, de-duplicates and sorts the read locations, and then estimates the
% genomic density signal using a variety of fixed- and adaptive-bandwidth kernel
% density estimation (KDE) techniques.
% 
% NOTES: 
% --------
%   (1) To save computational time during multiple function calls for processing
%       reads from the same dataset, the software stores the processed reads as .MAT files
%       in a separate subdirectory inside the current directory when a BED file is
%       processed for the first time. Two MAT files are created for each chromosome
%       corresponding to the positive- and the negative-strand reads, respectively. For
%       subsequent function calls involving the same dataset, the software automatically
%       detects the existence of processed reads and uses these directly from the MAT
%       files.
% 
%   (2) After processing a BED file, the software prints a report in the command
%       window with the details of the number of unique reads found in each
%       chromosome and the total number of reads processed.
% 
% All subfunctions necessary to use this code are contained in this file. The basic
% call to the main function is as follows:
% 
% MainKDE(fullbedname, fraglen, chrname, rangestart, rangeend, kdetype, kerneltype, kernelparam)
% 
% INPUTs:
% ----------
% The input parameters are defined as:
% 
% fullbedname -> the full name, including the extension, of a BED file containing the
%                read locations for the dataset of interest. If this file resides in
%                the same directory as this code, then there is no need to include
%                its path. Otherwise, the full path to the file has to be included.
% 
% fraglen     -> the mean fragment length, either estimated using a suitable method,
%                or whatever value the user desires to input. An accurate estimate of
%                the mean fragment length can be obtained using our recently released
%                software "MaSC", downloadable from the same webpage from where this
%                code was obtained. For details, see 
%                "MaSC: Mappability-Sensitive Cross-Correlation to Estimate Fragment Length for 
%                Short Read Sequencing Data", Ramachandran P, Palidwor GA, Porter CJ,
%                Perkins TJ, Bioinformatics, 2013.
% 
% chrname     -> the name of the chromosome, specified as, for e.g., 'chr1', followed
%                by two other parameters, 'rangestart' and 'rangeend', specifying the
%                chromosomal start and end positions defining the range of interest
%                for computing the density estimate. If it is desired to compute the
%                estimate for an entire chromosome, then specify empty matrices in
%                place of 'rangestart' and 'rangeend', like this [].
% 
% kdetype     -> the kde type, specify either 'fixed' or 'adaptive'
% 
% kerneltype  -> the kernel type, specify one of 'gauss', 'triangle', or 'square'
% 
% kernelparam -> the kernel parameter; If kerneltype is 'fixed', then this parameter
%                corresponds to the bandwidth or the standard deviation (sigma) of
%                the kernel measured in terms of number of base pairs. If kerneltype
%                is 'adaptive', then this parameter corresponds to the k in the
%                k-nearest neighbor rule, which is used to determine the bandwidth
%                for each read location (the adaptive-bandwidth rule). According to
%                this rule, the bandwidth associated with each read location is equal
%                to the greatest of the distance to the kth nearest location and the
%                distance to the 1-nearest neighbor on either side. 
% 
% [OPTIONAL] 'ploton' -> This is an optional input that can be specified after
%                        kernelparam. If 'ploton' is specified, then the computed
%                        kernel density is plotted and stored as a MATLAB fig file in
%                        the current directory. However, if nothing is specified
%                        after kernelparam, then the software defaults to NOT
%                        plotting and storing a figure.
% 
% OUTPUTs:
% ----------
% .wig file -> The computed density values are output in the form of a "fixedStep"
%              .wig (wiggle) file 
%              (see http://genome.ucsc.edu/goldenPath/help/wiggle.html)
% 
% plot      -> Optionally, if 'ploton' as described above is specified, then a MATLAB
%              .fig file is created in the current directory.
%  
% WARNINGs: 
% ----------
%           (1) It is advisable to choose 'ploton' only when the chosen range is of
%           reasonable size, say, less than 50000 bp or so. If the range is larger or
%           stretches for an entire chromosome, then plotting and storing the figure
%           can get very time-consuming and cumbersome, and would require a large
%           amount of hard-drive space!
% 
%           (2) If the density for multiple chromosomes or multiple ranges need to be
%           computed, the function must be called multiple times, once for each
%           chromosome or range. This might, potentially, change in the future!
% 
%           (3) The user must be aware that depending on the size of the dataset and
%           the specified range, the computations may take anywhere between a few
%           minutes and a few hours. During our benchmarking tests, for human chr22
%           of a ChIP-Seq dataset containing around 100K reads, with adaptive
%           Gaussian kernel and 5-nearest neighbor bandwidths, it took 0.38 seconds
%           to compute bandwidths and 47.5 seconds to compute the density on a single
%           core of a SunFire x2250 computer running Linux with 32Gb RAM.  A more
%           efficient implementation may be developed in the future.
% 
% EXAMPLES:
% -------------
% (1) MainKDE('abc.bed', 100, 'chr1', 4000000, 4001000, 'fixed', 'gauss', 30, 'ploton')
%     The above call will compute the density for a specified range of chr1, and
%     since the kdetype is 'fixed', the kernelparam value of 30 will be taken as the
%     bandwidth. And since 'ploton' is specified, a plot of the computed density will
%     be stored as a MATLAB .fig file in addition to the usual .wig file.
% 
% (2) MainKDE('abc.bed', 100, 'chr1', [], [], 'adaptive', 'triangle', 5)
%     The above call will compute the density for the entire chr1 (NOTE: this might
%     take a while and the output file would be large!), and since the kdetype is
%     'adaptive', the kernelparam value of 5 would be taken as the k for the
%     k-nearest neighbor rule to determine bandwidths. 
% 
% **************************************************************************************

function MainKDE(fullbedname, fraglen, chrname, rangestart, rangeend, kdetype, kerneltype, kernelparam, varargin)

[fpath, fname, fext] = fileparts(fullbedname);

if size(varargin,2) == 0
    ploton = false;
elseif size(varargin,2) == 1 && strcmpi(varargin,'ploton')
    ploton = true;
end

disp('**********************************************************')
disp(['Processing file : ' fname])
if isempty(rangestart) && isempty(rangeend)
    disp(['Whole ' chrname])
elseif ~isempty(rangestart) && ~isempty(rangeend)
    disp([chrname ' | Range: ' num2str(rangestart) ' - ' num2str(rangeend)])
end
disp('**********************************************************')

% Get the reads in sorted order
if isdir(fname)
    % If the directory exists, then read from the MAT files
    disp('Processed reads already exist, so using them ...')
    [CurSrtdPlusRdStrts, CurSrtdPlusRdEnds, CurSrtdMinusRdStrts, CurSrtdMinusRdEnds] = getchrreads(fname, chrname, rangestart, rangeend);
else
    % Else, first process the raw BED file, create the MAT files, and then read from the MAT files
    disp('No Processed reads found! So processing the BED file now ...')
    tic
    ReadBED(fullbedname);
    toc
    disp('Getting the required reads from the stored MAT files ...')
    [CurSrtdPlusRdStrts, CurSrtdPlusRdEnds, CurSrtdMinusRdStrts, CurSrtdMinusRdEnds] = getchrreads(fname, chrname, rangestart, rangeend);
end

% Shifting the Plus-Starts and Minus-Ends towards each other and combining reads
CmBndRds = sort([CurSrtdPlusRdStrts+round(fraglen/2) CurSrtdMinusRdEnds-round(fraglen/2)]);
if ploton && ~isempty(rangestart) && ~isempty(rangeend)
    PlotRds = CmBndRds( find(CmBndRds>=rangestart,1,'first') : find(CmBndRds<=rangeend,1,'last') );
elseif ploton && isempty(rangestart) && isempty(rangeend)
    PlotRds = CmBndRds(1):CmBndRds(end);
end
if ~isempty(rangestart) && ~isempty(rangeend)
    rangevals = rangestart:rangeend;
else
    rangevals = CmBndRds(1):CmBndRds(end);
end

if strcmpi(kdetype,'fixed')
    disp('Computing the chosen kernel density estimate ...')
    tic
    switch lower(kerneltype)
        case 'gauss'
            dNsiTy = kdeFxdBW_Gauss2(CmBndRds, kernelparam, rangevals);
        case 'square'
            dNsiTy = kdeFxdBW_Sq(CmBndRds, kernelparam, rangevals);
        case 'triangle'
            dNsiTy = kdeFxdBW_Tri(CmBndRds, kernelparam, rangevals);
        otherwise
            disp('ERROR! Unexpected kernel type. No density created. Enter one of gauss, triangle, or square within single quotes.'), beep
    end
    toc
elseif strcmpi(kdetype,'adaptive')
    tic
    disp('Computing the bandwidth corresponding to each read ...')
    LT = length(CmBndRds);
    BWs = zeros(1,LT);
    for jis = 1:LT
        CnDidteBWs = sort(abs(CmBndRds(jis) - CmBndRds(max(1,jis-kernelparam) : min(jis+kernelparam,LT))), 'ascend');
        OneFarBW = max(CmBndRds(jis)-CmBndRds(max(1,jis-1)), CmBndRds(min(LT,jis+1))-CmBndRds(jis));
        BWs(jis) = max(OneFarBW, CnDidteBWs(kernelparam+1));
    end
    toc
    disp('Computing the chosen kernel density estimate ...')
    tic
    switch lower(kerneltype)
        case 'gauss'
            dNsiTy = kdeVarBW_Gauss2(CmBndRds, BWs, rangevals);
        case 'square'
            dNsiTy = kdeVarBW_Sq(CmBndRds, BWs, rangevals);
        case 'triangle'
            dNsiTy = kdeVarBW_Tri(CmBndRds, BWs, rangevals);
        otherwise
            disp('ERROR! Unexpected kernel type. No density created. Enter one of gauss, triangle, or square within single quotes.'), beep
    end
    toc
else
    disp('ERROR! Unexpected kde type. No density created. Enter either fixed or adaptive within single quotes.'), beep
end

% WRITE to the output file in 'wiggle' format
disp('Creating output file ...')
tic
if ~isempty(rangestart) && ~isempty(rangeend)
    % Create the header line
    fid = fopen([fname '_FL_' num2str(fraglen) '_' kdetype '_' kerneltype '_param_' num2str(kernelparam) '_' ...
        chrname '_' num2str(rangestart) '_' num2str(rangeend) '.wig'],'w');
    fprintf(fid, 'fixedStep chrom=%s start=%d step=1\n', chrname, rangestart);
    fclose(fid);
    % End header line creation
    dlmwrite([fname '_FL_' num2str(fraglen) '_' kdetype '_' kerneltype '_param_' num2str(kernelparam) '_' ...
        chrname '_' num2str(rangestart) '_' num2str(rangeend) '.wig'], dNsiTy(:),'newline','pc','precision','%20.19f','-append');
else
    % Create the header line
    fid = fopen([fname '_FL_' num2str(fraglen) '_' kdetype '_' kerneltype '_param_' num2str(kernelparam) '_' ...
        chrname '.wig'],'w');
    fprintf(fid, 'fixedStep chrom=%s start=%d step=1\n', chrname, CmBndRds(1));
    fclose(fid);
    % End header line creation
    dlmwrite([fname '_FL_' num2str(fraglen) '_' kdetype '_' kerneltype '_param_' num2str(kernelparam) '_' ...
        chrname '.wig'], dNsiTy(:),'newline','pc','precision','%20.19f','-append');
end
toc

% Plotting, if chosen
disp('Creating plot ...')
tic
if ploton
    h1 = figure;
    plot(rangevals, dNsiTy),
    hold on
    YAxLims = get(gca, 'ylim');
    set(gca,'ylim',[YAxLims(1) 1.1*max(dNsiTy)]),
    Y1 = YAxLims(1)+0.05*diff(YAxLims);
    Y2 = YAxLims(1)+0.08*diff(YAxLims);
    for jjii = 1:length(PlotRds)
        line([PlotRds(jjii) PlotRds(jjii)], [Y1 Y2], 'color', 'r');
    end
    xlabel('Chromosomal locations')
    ylabel('Probability density')
    saveas(h1, [fname '_FL_' num2str(fraglen) '_' kdetype '_' kerneltype '_param_' num2str(kernelparam) '_' ...
        chrname '_' num2str(rangestart) '_' num2str(rangeend) '.fig'], 'fig');
    close all
end
toc
disp('All jobs complete!')


% Read a BED file
% separate into positive- and negative-strand reads, de-duplicate the reads, sort,
% and store for future use in a separate folder in the current directory as MAT files
function ReadBED(fullbedname)

[fpath, fname, fext] = fileparts(fullbedname);
fid = fopen(fullbedname);

disp('Loading the BED file ...')
% tTemp = fgetl(fid);  % Uncomment to reject the first line if this line has header info
MainCell = textscan(fid, '%s %f %f %*s %*s %s %*[^\n]', 'commentstyle', 't');  % Read only 4 columns and ignore the rest

status = fclose(fid);
if status ~= 0
    disp('The BED file did not close properly!')
end

% Just getting the total number of reads read
szMain = size(MainCell{2});
TotNoRds = szMain(1);
disp(['Total number of reads found in the BED file : ' num2str(TotNoRds)]);

% Obtaining the list of chromosomes present using "unique"
disp('Obtaining the list of chromosomes present ...')
[chrLst,~,chrIdx] = unique(MainCell{1});

disp('Separating reads, removing duplicates, and sorting ...')
AccumRds = struct;
for jis = 1:length(chrLst)
    % Plus
    Indx = (chrIdx == jis) & strcmpi(MainCell{4}, '+');
    [uPlusRdStrts,uIndx,~] = unique(MainCell{2}(Indx));
    noNuPlusRdEnds = MainCell{3}(Indx);
    uPlusRdEnds = noNuPlusRdEnds(uIndx);
    AccumRds.([chrLst{jis} '_Plus']).RdStrts = uPlusRdStrts;
    AccumRds.([chrLst{jis} '_Plus']).RdEnds  = uPlusRdEnds;
    %     Use the two lines below if the first colmn in the BED file has only nos. like "1"
    %     instead of "chr1"
    %     AccumRds.(['chr' chrLst{jis} '_Plus']).RdStrts = uPlusRdStrts;
    %     AccumRds.(['chr' chrLst{jis} '_Plus']).RdEnds  = uPlusRdEnds;
    
    % Minus
    Indx = (chrIdx == jis) & strcmpi(MainCell{4}, '-');
    [uMinusRdStrts,uIndx,~] = unique(MainCell{2}(Indx));
    noNuMinusRdEnds = MainCell{3}(Indx);
    uMinusRdEnds = noNuMinusRdEnds(uIndx);
    AccumRds.([chrLst{jis} '_Minus']).RdStrts = uMinusRdStrts;
    AccumRds.([chrLst{jis} '_Minus']).RdEnds  = uMinusRdEnds;
    %     Use the two lines below if the first colmn in BED file has only nos. like "1"
    %     instead of "chr1"
    %     AccumRds.(['chr' chrLst{jis} '_Minus']).RdStrts = uMinusRdStrts;
    %     AccumRds.(['chr' chrLst{jis} '_Minus']).RdEnds  = uMinusRdEnds;
    fprintf('%s', [chrLst{jis} ' '])
end
disp(' ')
disp('Deleting existing MAT files, if any, before saving new files ...')
if isdir(fname)
    rmdir(fname,'s');
end
mkdir(fname)
wkgtmpdir = cd(fname);
disp('Saving to MAT files ...')
AccumRds = orderfields(AccumRds);
fldList = fieldnames(AccumRds);
savedCnt = 0;
for jisw = 1:length(fldList)
    CurfldNam = fldList{jisw};
    if ~isempty(AccumRds.(CurfldNam).RdStrts)
        RdStrts = AccumRds.(CurfldNam).RdStrts(:)';
        RdEnds  = AccumRds.(CurfldNam).RdEnds(:)';
        save([fname '_' CurfldNam '_UniqSrtd_MatFile.mat'], 'RdStrts', 'RdEnds');
        savedCnt = savedCnt + 1;
    else
        disp(['No reads found for ' CurfldNam ' !!!'])
    end
end
cd(wkgtmpdir);

disp(' '), disp('*************** BED File Processing Report **************'), disp(' ')
disp(['Total number of files saved: ' num2str(savedCnt)]), disp(' ')
for jisw = 1:length(fldList)
    CurfldNam = fldList{jisw};
    Report(jisw,:) = {CurfldNam length(AccumRds.(CurfldNam).RdStrts)};
end
% Adding header line for report
Report = [{'chr Name'} {'No. of unique Reads found'}; {'*********'} {'*********************************'}; Report];
disp(Report)
disp(['Total number of Unique Reads found: ' num2str(sum(cell2mat(Report(3:end,2))))]),


% Get the reads for a particular chromosome (either whole or within a specified range) in sorted order
function [CurSrtdPlusRdStrts, CurSrtdPlusRdEnds, CurSrtdMinusRdStrts, CurSrtdMinusRdEnds] = getchrreads(fname, chrname, rangestart, rangeend)

oldfldr = cd(fname);

filLst = dir('*.mat');

% Locate the correct files, both Plus and Minus
foundPFile = false;  foundMFile = false;
for jis1 = 1:length(filLst)
    if ~isempty(strfind(filLst(jis1).name, [chrname '_Plus']))
        PlusFLNm = filLst(jis1).name;  foundPFile = true;
    end
    if ~isempty(strfind(filLst(jis1).name, [chrname '_Minus']))
        MinusFLNm = filLst(jis1).name;  foundMFile = true;
    end
end
if ~(foundPFile && foundMFile)
    disp('ERROR! At least one matching file, (Plus or Minus), NOT found for the selected chromosome!'), beep
end

% Load the two files and get the reads
PlusStruct  = load(PlusFLNm);
MinusStruct = load(MinusFLNm);

if isempty(rangestart) && isempty(rangeend)
    % Plus Reads
    PlusRdStrts = PlusStruct.RdStrts;
    PlusRdEnds  = PlusStruct.RdEnds;
    % Minus Reads
    MinusRdStrts = MinusStruct.RdStrts;
    MinusRdEnds  = MinusStruct.RdEnds;
elseif ~isempty(rangestart) && ~isempty(rangeend)
    % Plus Reads within range
    Plusind = PlusStruct.RdStrts >= rangestart & PlusStruct.RdStrts <= rangeend & PlusStruct.RdEnds >= rangestart & PlusStruct.RdEnds <= rangeend;
    PlusRdStrts = PlusStruct.RdStrts(Plusind);
    PlusRdEnds =  PlusStruct.RdEnds(Plusind);
    % Minus Reads within range
    Minusind = MinusStruct.RdStrts >= rangestart & MinusStruct.RdStrts <= rangeend & MinusStruct.RdEnds >= rangestart & MinusStruct.RdEnds <= rangeend;
    MinusRdStrts = MinusStruct.RdStrts(Minusind);
    MinusRdEnds =  MinusStruct.RdEnds(Minusind);
end

% Just to be sure that the reads are sorted
if ~issorted(PlusRdStrts)
    [CurSrtdPlusRdStrts, Pindx]   = sort(PlusRdStrts,  'ascend');
    CurSrtdPlusRdEnds  = PlusRdEnds(Pindx);
else
    CurSrtdPlusRdStrts = PlusRdStrts;
    CurSrtdPlusRdEnds  = PlusRdEnds;
end
if ~issorted(MinusRdStrts)
    [CurSrtdMinusRdStrts, Mindx]  = sort(MinusRdStrts,  'ascend');
    CurSrtdMinusRdEnds = MinusRdEnds(Mindx);
else
    CurSrtdMinusRdStrts = MinusRdStrts;
    CurSrtdMinusRdEnds = MinusRdEnds;
end
cd(oldfldr);

% This function carries out fixed BW KDE using Gaussian kernel
function dNsiTy = kdeFxdBW_Gauss2(srtdReads, BW, CurxMshRange)
% NOTE: Here, BW = sigma
n = length(srtdReads);
N = length(CurxMshRange);
dNsiTy = zeros(1,N);
denn = n;
halfshft = 5*BW; % We decide to truncate the kernel like this
sigma = BW;
for k = 1:n
    ReltvReadStrt = srtdReads(k)-CurxMshRange(1)+1;
    WnDw = max(1,ReltvReadStrt-halfshft) : min(N,ReltvReadStrt+halfshft);
    kern = exp(-0.5 * ((CurxMshRange(WnDw) - srtdReads(k))./sigma).^2) ./ sigma;
    kern = kern/sum(kern);
    dNsiTy(WnDw) = dNsiTy(WnDw) + kern;
end
% Normalizing
dNsiTy = dNsiTy/denn;

% This function carries out fixed BW KDE using square kernel
function dNsiTy  = kdeFxdBW_Sq(srtdReads, BW, CurxMshRange)
% NOTE: Here, BW = sigma
n       = length(srtdReads);
N       = length(CurxMshRange);
dNsiTy  = zeros(1,N);
% Unit Step kernel function
halfshft = floor(sqrt(3)*BW);
% The above sqrt(3) serves to get the equivalent std dev. as the Gaussian kernel
for k = 1:n
    ReltvReadStrt = srtdReads(k)-CurxMshRange(1)+1;
    WnDw = max(1,ReltvReadStrt-halfshft) : min(N,ReltvReadStrt+halfshft);
    dNsiTy(WnDw) = dNsiTy(WnDw) + 1/length(WnDw);
end
% Normalizing by the total number of reads
dNsiTy = dNsiTy/n;

% This function carries out Fxd BW KDE using Triangular kernel
function dNsiTy  = kdeFxdBW_Tri(srtdReads, BW, CurxMshRange)
% NOTE: Here, BW = sigma
n       = length(srtdReads);
N       = length(CurxMshRange);
dNsiTy  = zeros(1,N);
%  Triangular kernel function
halfshft = floor(sqrt(6)*BW);
% The above sqrt(6) is to get the equivalent std dev. as the Gaussian kernel
for k = 1:n
    ReltvReadStrt = srtdReads(k)-CurxMshRange(1)+1;
    if ReltvReadStrt == 1
        WnDw2 = ReltvReadStrt+1 : min(N,ReltvReadStrt+halfshft);
        c = CurxMshRange(ReltvReadStrt);
        b = CurxMshRange(WnDw2(end));
        kern = (2/(length(WnDw2)))*(b-CurxMshRange(WnDw2))/(b-c);
        kern = kern/sum(kern);
        dNsiTy(WnDw2) = dNsiTy(WnDw2) + kern;
    elseif ReltvReadStrt == N
        WnDw1 = max(1,ReltvReadStrt-halfshft) : ReltvReadStrt;
        a = CurxMshRange(WnDw1(1));
        c = CurxMshRange(ReltvReadStrt);
        kern = (2/(length(WnDw1)))*(CurxMshRange(WnDw1)-a)/(c-a);
        kern = kern/sum(kern);
        dNsiTy(WnDw1) = dNsiTy(WnDw1) + kern;
    else
        WnDw1 = max(1,ReltvReadStrt-halfshft) : ReltvReadStrt;
        WnDw2 = ReltvReadStrt+1 : min(N,ReltvReadStrt+halfshft);
        a = CurxMshRange(WnDw1(1));
        c = CurxMshRange(ReltvReadStrt);
        b = CurxMshRange(WnDw2(end));
        kern = (2/(length(WnDw1)+length(WnDw2)))*[(CurxMshRange(WnDw1)-a)/(c-a) (b-CurxMshRange(WnDw2))/(b-c)];
        kern = kern/sum(kern);
        dNsiTy([WnDw1 WnDw2]) = dNsiTy([WnDw1 WnDw2]) + kern;
    end
end
% Normalizing by the total number of reads
dNsiTy = dNsiTy/n;


% This function carries out Var BW KDE using Gaussian kernel
function dNsiTy = kdeVarBW_Gauss2(srtdReads, BWs, CurxMshRange)
% NOTE: Here, BW = sigma
n = length(srtdReads);
N = length(CurxMshRange);
dNsiTy = zeros(1,N);
LeftPos = zeros(1,n); LeftDrop = zeros(1,n);
RightPos = zeros(1,n); RightDrop = zeros(1,n);
denn = n;
for k = 1:n
    halfshft = 5*BWs(k); % We decide to truncate the kernel like this
    sigma = BWs(k);
    ReltvReadStrt = srtdReads(k)-CurxMshRange(1)+1;
    WnDw = max(1,ReltvReadStrt-halfshft) : min(N,ReltvReadStrt+halfshft);
    kern = exp(-0.5 * ((CurxMshRange(WnDw) - srtdReads(k))./sigma).^2) ./ sigma;
    kern = kern/sum(kern);
    dNsiTy(WnDw) = dNsiTy(WnDw) + kern;
end
% Normalizing
dNsiTy = dNsiTy/denn;


% This function carries out Var BW KDE using Square kernel
function dNsiTy  = kdeVarBW_Sq(srtdReads, BWs, CurxMshRange)
% NOTE: Here, BW = sigma
n       = length(srtdReads);
N       = length(CurxMshRange);
dNsiTy  = zeros(1,N);
for k = 1:n
    halfshft = floor(sqrt(3)*BWs(k));
    % The above sqrt(3) is to get the equivalent std dev. as the Gaussian kernel
    ReltvReadStrt = srtdReads(k)-CurxMshRange(1)+1;
    WnDw = max(1,ReltvReadStrt-halfshft) : min(N,ReltvReadStrt+halfshft);
    dNsiTy(WnDw) = dNsiTy(WnDw) + 1/length(WnDw);
end
% Normalizing by the total number of reads
dNsiTy = dNsiTy/n;


% This function carries out Var BW KDE using Triangular kernel
function dNsiTy = kdeVarBW_Tri(srtdReads, BWs, CurxMshRange)
% NOTE: Here, BW = sigma
n       = length(srtdReads);
N       = length(CurxMshRange);
dNsiTy  = zeros(1,N);
for k = 1:n
    halfshft = floor(sqrt(6)*BWs(k));
    % The above sqrt(6) is to get the equivalent std dev. as the Gaussian kernel
    ReltvReadStrt = srtdReads(k)-CurxMshRange(1)+1;
    if ReltvReadStrt == 1
        WnDw2 = ReltvReadStrt+1 : min(N,ReltvReadStrt+halfshft);
        c = CurxMshRange(ReltvReadStrt);
        b = CurxMshRange(WnDw2(end));
        kern = (2/(length(WnDw2)))*(b-CurxMshRange(WnDw2))/(b-c);
        kern = kern/sum(kern);
        dNsiTy(WnDw2) = dNsiTy(WnDw2) + kern;
    elseif ReltvReadStrt == N
        WnDw1 = max(1,ReltvReadStrt-halfshft) : ReltvReadStrt;
        a = CurxMshRange(WnDw1(1));
        c = CurxMshRange(ReltvReadStrt);
        kern = (2/(length(WnDw1)))*(CurxMshRange(WnDw1)-a)/(c-a);
        kern = kern/sum(kern);
        dNsiTy(WnDw1) = dNsiTy(WnDw1) + kern;
    else
        WnDw1 = max(1,ReltvReadStrt-halfshft) : ReltvReadStrt;
        WnDw2 = ReltvReadStrt+1 : min(N,ReltvReadStrt+halfshft);
        a = CurxMshRange(WnDw1(1));
        c = CurxMshRange(ReltvReadStrt);
        b = CurxMshRange(WnDw2(end));
        kern = (2/(length(WnDw1)+length(WnDw2)))*[(CurxMshRange(WnDw1)-a)/(c-a) (b-CurxMshRange(WnDw2))/(b-c)];
        kern = kern/sum(kern);
        dNsiTy([WnDw1 WnDw2]) = dNsiTy([WnDw1 WnDw2]) + kern;
    end
end
% Normalizing by the total number of reads
dNsiTy = dNsiTy/n;

