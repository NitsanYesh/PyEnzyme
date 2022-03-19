# Generated by PySCeS 1.0.1 (2022-03-17 22:47)
 
# Keywords
Description: 3IZNOK_TEST
Modelname: 3IZNOK_TEST
Output_In_Conc: True
Species_In_Conc: True
 
# GlobalUnitDefinitions
UnitSubstance: mole, 1.0, 0, 1
UnitVolume: litre, 1.0, 0, 1
UnitTime: second, 1.0, 0, 1
UnitLength: metre, 1.0, 0, 1
UnitArea: metre, 1.0, 0, 2
 
FIX: p0 
 
# Compartments
Compartment: v0, nan, 3 
 
# Reactions
r0@v0:
    s0 + s1 > s2 + s3
    r0_k_cat*p0*s1/(r0_k_m+s1)
# r0 has modifier(s): p0  
 
# Fixed species
p0@v0 = nan
 
# Variable species
s0@v0 = nan
s1@v0 = nan
s2@v0 = nan
s3@v0 = nan
 
# Parameters
r0_k_cat = 0.015
r0_k_m = 0.01
 
