# Distributed Approach to Satellite Direct-to-Cell Connectivity in 6G Non-Terrestrial Networks
This repository provides the code generated for the following article:
- Title: Distributed Approach to Satellite Direct-to-Cell Connectivity in 6G Non-Terrestrial Networks
- Authors: Diego Tuzi, Estephania Flores Aguilar, Thomas Delamotte, Gunes Karabulut-Kurt, and Andreas Knopp
- Link: https://ieeexplore.ieee.org/document/10355084

## Description
The code has been written using MATLAB R2023b and it is basically subdivided into two files:
- main.m
- figures.m

The function elsaGeometry.m implements the enhanced logarithmic spiral array (ELSA) geometry described in [Satellite Swarm-Based Antenna Arrays for 6G Direct-to-Cell Connectivity](https://ieeexplore.ieee.org/document/10068542).

## Main usage
Option 1:
1. Customize your parameters in the first lines of the main.m script
2. Run the main.m script with MATLAB. Results are saved into the "data" folder
3. Run the figures.m script with MATLAB loading the results of the previous point

Option 2:
1. Run the figures.m script with MATLAB loading the default results provided in data/results_article.mat

## Contacts
Authors: 
- [Diego Tuzi](https://www.linkedin.com/in/diegotuzi/) (corresponding author), Thomas Delamotte, and Andreas Knopp are with the University of the Bundeswehr Munich, Germany
- Estephania Flores Aguilar and Gunes Karabulut-Kurt are with Polytechnique Montréal, Canada

Email: <paper.sp@unibw.de>
