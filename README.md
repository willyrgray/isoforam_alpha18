ISOFORAM alpha-18 - species specific O-18 fractionation relationships for planktic and benthic foraminifera

solves for d18O_calcite/temperature/d18O_sw given the other two as inputs using foram species-specific fractionation factors from Daeron + Gray 2023

Revisiting 18-O and clumped isotopes in Planktic and Benthic foraminifera https://doi.org/10.1029/2023PA004660 

if flag is 1, returns d18Oc on VPDB scale given var1 = T in celsius and var2 is d18Osw on VSMOW

if flag is 2, returns T in celcius given var1 = d18Oc on VPDB scale and var2 is d18Osw on VSMOW

if flag is 3, returns d18Osw on VSMOW scale given var1 = T in celsius and var2 is d18Oc on VPDB scale

VSMOW VPDB conversion is IUAPC (Brand 2014 / Kim 2015)

isoforam_alpha18() returns central value plus calibration uncertainites (SE and SD)

isoforam_a18() returns central value plus calibration uncertainites (SE and SD)

for 'species' input the number corresponding to the species/genus/group listed below
i.e. for 'All planktics' input species=33 
[1] "Globigerina bulloides"            
[2] "Globigerinoides ruber white"      
[3] "Globorotalia inflata"             
[4] "Globorotalia truncatulinoides (d)"
[5] "Globorotalia truncatulinoides (s)"
[6] "Neogloboquadrina pachyderma"      
[7] "Orbulina universa"                
[8] "Trilobatus sacculifer"            
[9] "Bulimina aculeata"                
[10] "Bulimina marginata"               
[11] "Cibicidoides pachyderma"          
[12] "Cibicidoides wullerstorfi"        
[13] "Heoglundina elegans"              
[14] "Planulina ariminensis"            
[15] "Planulina foveolata"              
[16] "Rosalina vilardeboana"            
[17] "Uvigerina curticosta"             
[18] "Uvigerina flintii"                
[19] "Uvigerina peregrina"              
[20] "Bulimina"                         
[21] "Cibicidoides"                     
[22] "Globigerina"                      
[23] "Globigerinoides"                  
[24] "Globorotalia"                     
[25] "Hoeglundina"                      
[26] "Neogloboquadrina"                 
[27] "Orbulina"                         
[28] "Planulina"                        
[29] "Rosalina"                         
[30] "Trilobatus"                       
[31] "Uvigerina"                        
[32] "Cibicidoides + Planulina"         
[33] "All planktics"
