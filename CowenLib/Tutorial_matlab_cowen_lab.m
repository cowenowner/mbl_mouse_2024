%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a general matlab tutorial of the fundamentals needed for any
% analysis of neural data. More specific tutorials for analyzing specific
% types of neural data can be found elsewhere, but without the foundations,
% you will have trouble understanding, debugging, and expanding your
% analyses.
% 
% This tutorial is 'quest based' in that it poses objectives without
% answers and then the matlab adventurer is charged with completing the
% quest. Book learning for matlab only gets you so far. You need to type
% stuff, err, and try and err again. Hints are given in parentheses for
% each quest! Feel free to use the matlab help (doc <function>) or the web
% or Chat-GPT or friends. BUT, do your own work or you won't learn.
%
% The first step of this tutorial is to create a copy of this file so that
% you can edit it. Rename it so that your name is at the end. For example
% Tutorial_matlab_cowen_lab_ChuckieCheese.m
%
% This tutorial presumes you know how to start matlab and to create your
% own .m file (a text file that stores your code). It also assumes that you
% know what it means to have files/functions in your matlab path. If you do
% not know these things, do some background research before starting.
%
% Cowen 2024

% Q1: Create a histogram of 1000 random numbers with a mean of 10 and a
% standard deviation of 4. (randn, histogram)

% Q2: Take your list of 1000 random numbers - also called a vector - and 1)
% assign it to a variable. So that we're all consistent, use the variable
% name rnd_nums. 2) add 10 3) divide by 14, 4) multiply by 20, and 5) do an
% element by element multiplication of that list with itself. (hint: add a
% dot before the operator. An 'operator' is something like +,-,*,/).

% Q3: Make 2 different vectors of random numbers with mean 5 and std of 2.
% Assign 1 to the variable named V1 and another to the variable named V2.
% Plot a scatter plot (plot function with parameters such as '.' or 'o' or
% use the scatter function instead) of V1 on the x axis and V2 on the y
% axis.

% Q4: Matrices and extracting single columns and rows: 1) create a random
% matrix of 100 rows and 2000 columns using the rand function in matlab.
% Use the help. Assign that matrix to the variable M. 2) visualize that
% matrix using imagesc 3) plot a single row (row number 5) from that matrix
% using the plot function. 4) plot a single column - column 144 - using the
% plot function. 5) replot the matrix using imagesc, but this time
% transpose the matrix before plotting it. What is the difference?

% Q5: Create a vector of ascending numbers from 1 to 1000 but skip every 5
% numbers. (hint: 1:2:10)

% Q6: Create your first function. Matlab functions are different from
% scripts. Functions are just like math functions you know and love like
% sine() cosine(). You give the function some input (in paretheses) and it
% spits out an output: data = sine(0:.01:10*pi) will spit out a sine wave
% from radians ranging from 0 to 10:pi radians. Each function typically has
% it's own dedicated .m file. Create a function called my_dumb_function.m
% and save it somewhere. Make this function take in a single vector of
% numbers and spit out the sum cubed of this number. hint(sum, .^3). Then
% test your function by calling it on the command line. Be sure it is in
% your path.

% Now that you have these basics, use the help to look up these useful
% functions: isempty, readtable, writetable, addpath, genpath, dir, mkdir,
% cd, find. If you have some time, use these functions - write some code
% that uses them.
%
% even more important, look up some of the 'data types'. You will use these
% every day: double, single, cell array, structure, structure array, table.