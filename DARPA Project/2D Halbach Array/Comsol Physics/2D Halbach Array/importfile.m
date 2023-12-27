function [x, y, z, Bx, By, Bz] = importfile(filename)
%IMPORTFILE Import data from a text file
%  [X, Y, Z, BX, BY, BZ] = IMPORTFILE(FILENAME) reads data from text
%  file FILENAME for the default selection.  Returns the data as column
%  vectors.
%
%  [X, Y, Z, BX, BY, BZ] = IMPORTFILE(FILE, DATALINES) reads data for
%  the specified row interval(s) of text file FILENAME. Specify
%  DATALINES as a positive scalar integer or a N-by-2 array of positive
%  scalar integers for dis-contiguous row intervals.
%
%  Example:
%  [x, y, z, Bx, By, Bz] = importfile("D:\Box\GitHub\liquid_mirrors\Comsol_Simulations\B_0.5_0.5_0.5.txt", [10, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 17-Dec-2023 22:21:12

%% Input handling

% If dataLines is not specified, define defaults
dataLines = [10, Inf];

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 17);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["x", "y", "z", "Bx", "By", "Bz", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17"];
opts.SelectedVariableNames = ["x", "y", "z", "Bx", "By", "Bz"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["x", "y", "z", "Bx", "By", "Bz"], "ThousandsSeparator", ",");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
x = tbl.x;
y = tbl.y;
z = tbl.z;
Bx = tbl.Bx;
By = tbl.By;
Bz = tbl.Bz;
end