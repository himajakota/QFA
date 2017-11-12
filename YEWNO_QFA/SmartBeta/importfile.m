function tableout = importfile(workbookFile,sheetName,startRow,endRow)
%IMPORTFILE1 Import data from a spreadsheet
%   DATA = IMPORTFILE1(FILE) reads data from the first worksheet in the
%   Microsoft Excel spreadsheet file named FILE and returns the data as a
%   table.
%
%   DATA = IMPORTFILE1(FILE,SHEET) reads from the specified worksheet.
%
%   DATA = IMPORTFILE1(FILE,SHEET,STARTROW,ENDROW) reads from the specified
%   worksheet for the specified row interval(s). Specify STARTROW and
%   ENDROW as a pair of scalars or vectors of matching size for
%   dis-contiguous row intervals. To read to the end of the file specify an
%   ENDROW of inf.%
% Example:
%   priceData = importfile1('priceData.xlsx','Sheet1',2,254);
%
%   See also XLSREAD.

% Auto-generated by MATLAB on 2017/11/11 15:55:41

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    startRow = 2;
    endRow = 254;
end

%% Import the data, extracting spreadsheet dates in Excel serial date format
[~, ~, raw, dates] = xlsread(workbookFile, sheetName, sprintf('A%d:AC%d',startRow(1),endRow(1)),'' , @convertSpreadsheetExcelDates);
for block=2:length(startRow)
    [~, ~, tmpRawBlock,tmpDateNumBlock] = xlsread(workbookFile, sheetName, sprintf('A%d:AC%d',startRow(block),endRow(block)),'' , @convertSpreadsheetExcelDates);
    raw = [raw;tmpRawBlock]; %#ok<AGROW>
    dates = [dates;tmpDateNumBlock]; %#ok<AGROW>
end
raw = raw(:,[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]);
dates = dates(:,1);

%% Create output variable
I = cellfun(@(x) ischar(x), raw);
raw(I) = {NaN};
data = reshape([raw{:}],size(raw));

%% Create table
tableout = table;

%% Allocate imported array to column variable names
dates(~cellfun(@(x) isnumeric(x) || islogical(x), dates)) = {NaN};
tableout.Date = datetime([dates{:,1}].', 'ConvertFrom', 'Excel');
tableout.APPL = data(:,1);
tableout.MSFT = data(:,2);
tableout.PFE = data(:,3);
tableout.JNJ = data(:,4);
tableout.CAT = data(:,5);
tableout.INTC = data(:,6);
tableout.HD = data(:,7);
tableout.MCD = data(:,8);
tableout.UNH = data(:,9);
tableout.MMM = data(:,10);
tableout.GS = data(:,11);
tableout.V = data(:,12);
tableout.CVX = data(:,13);
tableout.TRV = data(:,14);
tableout.PG = data(:,15);
tableout.DIS = data(:,16);
tableout.WMT = data(:,17);
tableout.AXP = data(:,18);
tableout.JPM = data(:,19);
tableout.GE = data(:,20);
tableout.DWDP = data(:,21);
tableout.MRK = data(:,22);
tableout.NKE = data(:,23);
tableout.KO = data(:,24);
tableout.VZ = data(:,25);
tableout.CSCO = data(:,26);
tableout.BA = data(:,27);
tableout.SPDR = data(:,28);

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% tableout.Date=datenum(tableout.Date);

