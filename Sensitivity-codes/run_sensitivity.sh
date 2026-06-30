#!/bin/bash
#$ -M email@nd.edu (change this to receive notification about the simualtion start and end time)
#$ -m abe
#$ -q long

module load abaqus intel/24.2
abaqus cae noGUI=sensitivity_analysis.py
