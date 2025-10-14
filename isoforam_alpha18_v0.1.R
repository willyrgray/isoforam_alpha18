#ISOFORAM alpha-18 - species specific O-18 fractionation relationships for planktic and benthic foraminifera
#Daeron + Gray 2023
#Revisiting 18-O and clumped isotopes in Planktic and Benthic foraminifera https://doi.org/10.1029/2023PA004660 
#if flag is 1, returns d18Oc on VPDB scale given var1 = T in celsius and var2 is d18Osw on VSMOW
#if flag is 2, returns T in celcius given var1 = d18Oc on VPDB scale and var2 is d18Osw on VSMOW
#if flag is 3, returns d18Osw on VSMOW scale given var1 = T in celsius and var2 is d18Oc on VPDB scale
#VSMOW VPDB conversion is IUAPC (Brand 2014 / Kim 2015)

#for 'species' input the number corresponding to the species/genus/group listed below
#i.e. for 'All planktics' input species=33 
#[1] "Globigerina bulloides"            
#[2] "Globigerinoides ruber white"      
#[3] "Globorotalia inflata"             
#[4] "Globorotalia truncatulinoides (d)"
#[5] "Globorotalia truncatulinoides (s)"
#[6] "Neogloboquadrina pachyderma"      
#[7] "Orbulina universa"                
#[8] "Trilobatus sacculifer"            
#[9] "Bulimina aculeata"                
#[10] "Bulimina marginata"               
#[11] "Cibicidoides pachyderma"          
#[12] "Cibicidoides wullerstorfi"        
#[13] "Heoglundina elegans"              
#[14] "Planulina ariminensis"            
#[15] "Planulina foveolata"              
#[16] "Rosalina vilardeboana"            
#[17] "Uvigerina curticosta"             
#[18] "Uvigerina flintii"                
#[19] "Uvigerina peregrina"              
#[20] "Bulimina"                         
#[21] "Cibicidoides"                     
#[22] "Globigerina"                      
#[23] "Globigerinoides"                  
#[24] "Globorotalia"                     
#[25] "Hoeglundina"                      
#[26] "Neogloboquadrina"                 
#[27] "Orbulina"                         
#[28] "Planulina"                        
#[29] "Rosalina"                         
#[30] "Trilobatus"                       
#[31] "Uvigerina"                        
#[32] "Cibicidoides + Planulina"         
#[33] "All planktics"


#################
####start function 
isoforam_alpha18<- function(flag, species, var1, var2){
	#species=33
	if(is.null(species)){species == 33} #defaults to 'all planktics' if no species input
	
	#defined constants
	species_list<- c('Globigerina bulloides', 'Globigerinoides ruber white', 'Globorotalia inflata','Globorotalia truncatulinoides (d)', 'Globorotalia truncatulinoides (s)','Neogloboquadrina pachyderma','Orbulina universa','Trilobatus sacculifer','Bulimina aculeata','Bulimina marginata','Cibicidoides pachyderma','Cibicidoides wullerstorfi','Heoglundina elegans','Planulina ariminensis','Planulina foveolata','Rosalina vilardeboana','Uvigerina curticosta','Uvigerina flintii','Uvigerina peregrina','Bulimina','Cibicidoides','Globigerina','Globigerinoides','Globorotalia','Hoeglundina','Neogloboquadrina','Orbulina','Planulina','Rosalina','Trilobatus','Uvigerina','Cibicidoides + Planulina','All planktics') 
	
	B<- c(-32.47, -32.75, -32.05, -32.20, -32.10, -32.32, -31.95, -32.67, -31.87, -31.99, -32.22, -32.24, -31.06, -32.26, -32.32, -32.18, -31.40, -31.71, -31.73, -31.95, -32.23, -32.47, -32.75, -32.11, -31.06, -32.32, -31.95, -32.29, -32.18, -32.67, -31.73, -32.24, -32.49)
	
	B_SE<- c(0.031, 0.022, 0.063, 0.102, 0.049, 0.045, 0.041, 0.027, 0.020, 0.014, 0.025, 0.018, 0.026, 0.044, 0.056, 0.025, 0.037, 0.066, 0.039, 0.012, 0.013, 0.031, 0.022, 0.038, 0.026, 0.045, 0.041, 0.036, 0.025, 0.027, 0.025, 0.013, 0.020)
	
	B_SD<- c(0.21, 0.21, 0.23, 0.32, 0.23, 0.31, 0.18, 0.22, 0.13, 0.12, 0.13, 0.09, 0.27, 0.13, 0.18, 0.12, 0.1, 0.13, 0.2, 0.14, 0.11, 0.21, 0.21, 0.25, 0.27, 0.31, 0.18, 0.16, 0.12, 0.22, 0.21, 0.13, 0.35)
	
	###call constants
	A<- 18.03*10^3 #inorganic calcite value from kim o'neil
	
	Bspc<- B[species] #species specific offset from Kim'Oneil 97
	Bspc_se<- B_SE[species] #standard error of species specific offset 
	Bspc_sd<- B_SD[species] #standard deviation of species specific offset   
	#B<- -32.42 #inorganic calcite from kim o'neil

	if(flag == 1) {
		#returns d18Oc on VPDB scale given var1 = T in celsius and var2 is d18Osw on VSMOW
		#var1 = 1; var2 = 0
		alpha_ln1000<- A/(var1+273.15)+Bspc
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dc<- (var2+ (dc_dw - 30.92)/1.03092)
		#upr se
		alpha_ln1000<- A/(var1+273.15)+(Bspc+Bspc_se)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dc_se_upr<- (var2+ (dc_dw - 30.92)/1.03092)
		#lwr se
		alpha_ln1000<- A/(var1+273.15)+(Bspc-Bspc_se)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dc_se_lwr<- (var2+ (dc_dw - 30.92)/1.03092)
		#upr sd
		alpha_ln1000<- A/(var1+273.15)+(Bspc+Bspc_sd)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dc_sd_upr<- (var2+ (dc_dw - 30.92)/1.03092)
		#lwr sd
		alpha_ln1000<- A/(var1+273.15)+(Bspc-Bspc_sd)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dc_sd_lwr<- (var2+ (dc_dw - 30.92)/1.03092)
		#
		dc_se<- mean(c(abs(dc-dc_se_lwr), abs(dc_se_lwr-dc)))
		dc_sd<- mean(c(abs(dc-dc_sd_lwr), abs(dc_sd_lwr-dc)))
		return(data.frame(d18Oc=dc, d18Oc_se=dc_se, d18Oc_sd=dc_sd))
	} 

	else if(flag == 2) {
		#returns T in celcius given var1 = d18Oc on VPDB scale and var2 is d18Osw on VSMOW
		#var1 = 2.829303; var2 = 0					
		dc_dw_vsmow<-1.03092*(var1-var2)+30.92
		alpha<- dc_dw_vsmow/1000+1
		Tk<- A / ((log(alpha)*1000)-Bspc)				
		Tc<- Tk-273.15 
		#upr se
		dc_dw_vsmow<-1.03092*(var1-var2)+30.92
		alpha<- dc_dw_vsmow/1000+1
		Tk<- A / ((log(alpha)*1000)-(Bspc+Bspc_se))				
		Tc_se_upr<- Tk-273.15 
		#lwr se
		dc_dw_vsmow<-1.03092*(var1-var2)+30.92
		alpha<- dc_dw_vsmow/1000+1
		Tk<- A / ((log(alpha)*1000)-(Bspc-Bspc_se))				
		Tc_se_lwr<- Tk-273.15 
		#upr sd
		dc_dw_vsmow<-1.03092*(var1-var2)+30.92
		alpha<- dc_dw_vsmow/1000+1
		Tk<- A / ((log(alpha)*1000)-(Bspc+Bspc_sd))				
		Tc_sd_upr<- Tk-273.15 
		#lwr se
		dc_dw_vsmow<-1.03092*(var1-var2)+30.92
		alpha<- dc_dw_vsmow/1000+1
		Tk<- A / ((log(alpha)*1000)-(Bspc-Bspc_sd))				
		Tc_sd_lwr<- Tk-273.15 
		#
		Tc_se<- mean(c(abs(Tc-Tc_se_lwr), abs(Tc_se_lwr-Tc)))
		Tc_sd<- mean(c(abs(Tc-Tc_sd_lwr), abs(Tc_sd_lwr-Tc)))
		return(data.frame(t = Tc, t_se = Tc_se, t_sd = Tc_sd))
	}	
	
	else if(flag == 3) {
		#returns d18Osw on VSMOW scale given var1 = T in celsius and var2 is d18Oc on VPDB scale
		#var1 = 1; var2 = 2.829303
		alpha_ln1000<- A/(var1+273.15)+Bspc
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dw<- (var2-(dc_dw-30.92)/1.03092)
		#upr se
		alpha_ln1000<- A/(var1+273.15)+(Bspc+Bspc_se)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dw_se_upr<- (var2-(dc_dw-30.92)/1.03092)
		#lwr se
		alpha_ln1000<- A/(var1+273.15)+(Bspc-Bspc_se)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dw_se_lwr<- (var2-(dc_dw-30.92)/1.03092)
		#upr sd
		alpha_ln1000<- A/(var1+273.15)+(Bspc+Bspc_sd)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dw_sd_upr<- (var2-(dc_dw-30.92)/1.03092)
		#lwr sd
		alpha_ln1000<- A/(var1+273.15)+(Bspc-Bspc_sd)
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dw_sd_lwr<- (var2-(dc_dw-30.92)/1.03092)
		#
		dw_se<- mean(c(abs(dw-dw_se_lwr), abs(dw_se_lwr-dw)))
		dw_sd<- mean(c(abs(dw-dw_sd_lwr), abs(dw_sd_lwr-dw)))
		return(data.frame(d18Ow = dw, d18Ow_se = dw_se, d18Ow_sd = dw_sd))
	}
	#print(paste(species_list[species],'calibration')) 
}

###end function 
#################################

isoforam_alpha18(flag=1, species=32, var1=-2, var2=0)
# d18Oc   d18Oc_se  d18Oc_sd
# 3.810193 0.01304945 0.1304868

isoforam_alpha18(flag=2, species=32, var1=3.810193, var2=0)
#         t       t_se      t_sd
# -1.999999 0.05300074 0.5290767

#test function 
isoforam_alpha18(flag=2, species=32, var1=isoforam_alpha18(flag=1, species=32, var1=1, var2=0)$d18Oc, var2=0)
#t       t_se      t_sd
# 1 0.05417991 0.5408371 #it works!

isoforam_alpha18(flag=3, species=32, var1=-2, var2=3.810193)
#         d18Ow   d18Ow_se  d18Ow_sd
# -1.50902e-07 0.01304945 0.1304868


#test function 
isoforam_alpha18(flag=2, species=32, var1=isoforam_alpha18(flag=1, species=32, var1=1, var2=0)$d18Oc, var2=0)
#t       t_se      t_sd
# 1 0.05417991 0.5408371 #it works!




#############################################################################
#same function as above without returning errors to simplify output 
#################
####start function 
isoforam_a18<- function(flag, species, var1, var2){
	#species=33
	if(is.null(species)){species == 33} #defaults to 'all planktics' if no species input
	
	#defined constants
	species_list<- c('Globigerina bulloides', 'Globigerinoides ruber white', 'Globorotalia inflata','Globorotalia truncatulinoides (d)', 'Globorotalia truncatulinoides (s)','Neogloboquadrina pachyderma','Orbulina universa','Trilobatus sacculifer','Bulimina aculeata','Bulimina marginata','Cibicidoides pachyderma','Cibicidoides wullerstorfi','Heoglundina elegans','Planulina ariminensis','Planulina foveolata','Rosalina vilardeboana','Uvigerina curticosta','Uvigerina flintii','Uvigerina peregrina','Bulimina','Cibicidoides','Globigerina','Globigerinoides','Globorotalia','Hoeglundina','Neogloboquadrina','Orbulina','Planulina','Rosalina','Trilobatus','Uvigerina','Cibicidoides + Planulina','All planktics') 
	
	B<- c(-32.47, -32.75, -32.05, -32.20, -32.10, -32.32, -31.95, -32.67, -31.87, -31.99, -32.22, -32.24, -31.06, -32.26, -32.32, -32.18, -31.40, -31.71, -31.73, -31.95, -32.23, -32.47, -32.75, -32.11, -31.06, -32.32, -31.95, -32.29, -32.18, -32.67, -31.73, -32.24, -32.49)
	
	#B_SE<- c(0.031, 0.022, 0.063, 0.102, 0.049, 0.045, 0.041, 0.027, 0.020, 0.014, 0.025, 0.018, 0.026, 0.044, 0.056, 0.025, 0.037, 0.066, 0.039, 0.012, 0.013, 0.031, 0.022, 0.038, 0.026, 0.045, 0.041, 0.036, 0.025, 0.027, 0.025, 0.013, 0.020)
	
	#B_SD<- c(0.21, 0.21, 0.23, 0.32, 0.23, 0.31, 0.18, 0.22, 0.13, 0.12, 0.13, 0.09, 0.27, 0.13, 0.18, 0.12, 0.1, 0.13, 0.2, 0.14, 0.11, 0.21, 0.21, 0.25, 0.27, 0.31, 0.18, 0.16, 0.12, 0.22, 0.21, 0.13, 0.35)
	
	###call constants
	A<- 18.03*10^3 #inorganic calcite value from kim o'neil
	
	Bspc<- B[species] #species specific offset from Kim'Oneil 97
	#Bspc_se<- B_SE[species] #standard error of species specific offset 
	#Bspc_sd<- B_SD[species] #standard deviation of species specific offset   
	#B<- -32.42 #inorganic calcite from kim o'neil

	if(flag == 1) {
		#returns d18Oc on VPDB scale given var1 = T in celsius and var2 is d18Osw on VSMOW
		#var1 = 1; var2 = 0
		alpha_ln1000<- A/(var1+273.15)+Bspc
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dc<- (var2+ (dc_dw - 30.92)/1.03092)
		return(dc)
	} 
	
	else if(flag == 2) {
		#returns T in celcius given var1 = d18Oc on VPDB scale and var2 is d18Osw on VSMOW
		#var1 = 2.829303; var2 = 0					
		dc_dw_vsmow<-1.03092*(var1-var2)+30.92
		alpha<- dc_dw_vsmow/1000+1
		Tk<- A / ((log(alpha)*1000)-Bspc)				
		Tc<- Tk-273.15 
		return(Tc)
	}	
	
	else if(flag == 3) {
		#returns d18Osw on VSMOW scale given var1 = T in celsius and var2 is d18Oc on VPDB scale
		#var1 = 1; var2 = 2.829303
		alpha_ln1000<- A/(var1+273.15)+Bspc
		alpha<- exp(alpha_ln1000/1000)
		dc_dw<- (alpha-1)*1000
		dw<- (var2-(dc_dw-30.92)/1.03092)
		return(dw)
	}
	#print(paste(species_list[species],'calibration')) 
}

###end function 
#################################

#test function 
isoforam_a18(flag=2, species=32, var1=isoforam_a18(flag=1, species=32, var1=1, var2=0), var2=0)
#1
