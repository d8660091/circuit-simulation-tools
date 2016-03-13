%--------------------------------------------------------------------------
% aRunBlkArnoldi.m
% My MIMO Simulator using BLOCK ARNOLDI Algorithm
%--------------------------------------------------------------------------
close all;  clc;  clear all

addpath(genpath('SrvrFunc'));   %Call genpath inside of addpath to also add all subfolders
addpath(genpath('ArnoldiFunc'));   

dirFN.ORG = '../A1_ORG/!DO/MNA.mat';
dirFN.Q   = '!DO/Q.mat';
dirFN.H   = '!DO/H.mat';
dirFN.A   = '!DO/A.mat';

Nr_Block_Arnoldi_Moments = 15;

delete('!DO/Report_On_CW.*');
%diary on;
diary('!DO/Report_On_CW.txt');

load( dirFN.ORG, 'G', 'C', 'B', 'SF' );  %including: G,C,B,L, Nr.NrNodes, Nr.MNAsize, Nr.PortsNodes Nr.Ks

RunArnoldi( G, C, B, SF, dirFN, Nr_Block_Arnoldi_Moments )

xcput;
diary off;
fclose all;
%--------------------------------------------------------------------------


