# SimplusGT Coding Style

## Introduction

This is an instruction of the detailed coding style for Simplus Grid Tool (SimplusGT). Please read "[DeveloperManual.md](https://github.com/Future-Power-Networks/Simplus-Grid-Tool/blob/master/Documentations/DeveloperManual.md)" first before reading this file.

## Naming Convention

There are three sets of naming conventions:  
    1) microsoft, like PlotFigure  
    2) gnu-linux, like plot_figure  
    3) free style, like ss2tf  
Matlab uses free style for its common functions, and uses microsoft mainly in others. We also use microsoft for this toolbox, and use free style when needed. This rule is valid for names of folders, functions, variables.

## Folder

".m" Functions should be saved in the corresponding folders. All functions are saved in "+SimplusGT" folder, which is the root namespace. Generic functions (used for generic purpose such as mathmatical calculations, bode plot, etc) are saved in this root name space folder directly. Other advanced functions are saved in corresponding sub name space folders such as "+SimplusGT/+PowerFlow".

## Comment

There must be header comments at the beginning of each ".m" file, which briefly introduce this file and shows the original author(s) of this file. For example, in "+SimplusGT/abc2dq.m":

<pre>
% This function transforms a time-domain signal from the natural frame (abc
% frame) to synchronous frame (dq0 frame).
%
% Author(s): Yitong Li
</pre>

If others need to modify this file, please add the following statement, to indicate who modified the file and what was modified. For example:

<pre>
% Modified by: Yunjie Gu
% Modifications: ......
</pre>

Comments should also be added when writting the main codes of each ".m" file. Please use English to comment. Comments should be placed at the right of corresponding code normally.

## Indent

Please use TAB or 4 spaces to indent. For statements including "for", "if ... else ...", etc, indent is recommended. For example:

<pre>
for i = 1:n
    a = a + 1;
end
</pre>


