//ETHAN WHITMAN PACE ANALYSIS


//global db3275 "db3275"
global db3275 "danielbelsky"

cd "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Phenotype Data"
//Subject File
use "Subject&Pedigree_AllCohorts_v33.dta", clear 
//Ages & Dates
merge 1:1 dbgap_subject_id using "Age&Dates_AllCohorts_v33.dta", nogen
//Birth Year
merge 1:1 dbgap_subject_id using "BirthYear_AllCohorts_v33.dta", nogen
//Education
merge 1:1 dbgap_subject_id using "v33EducationData.dta", nogen
gen EA = 0 if education == "a:noHSdeg"
replace EA =1 if education == "b:HSde"
replace EA =2 if education == "c:somecoll"
replace EA =3 if education == "d:collgrad"
capture label drop EA 
label define EA 0 "<HS" 1 "HS" 2 "Some College" 3 "College Grad"
label values EA EA 
merge 1:1 dbgap_subject_id using "SES_AllCohorts_v33.dta", keepus(maxeducg maxeayrs) nogen 
//Smoking
merge 1:1 dbgap_subject_id using "SmokingStatus_OffspringGen3_v33.dta", nogen
//Survival
merge 1:1 dbgap_subject_id using "Survival_AllCohorts_v33.dta", nogen
//Dementia 
merge 1:1 dbgap_subject_id using "Dementia_AllCohorts_v33.dta", nogen 

//MRI data 
merge 1:1 dbgap_subject_id using "MRIOffspring2010_v33.dta", nogen 


cd "//Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2"

//DNAm clocks 
merge 1:1 dbgap_subject_id using "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/MolecularData/DNAm/Processed_Variables/Framingham_Offspring/Framingham_Offspring_DNAmVariables230628dwb.dta", nogen 

//APOE4 data
drop sample_id 
merge 1:1 dbgap_subject_id using "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/MolecularData/SHARe_imputed_phg000835.v5/APOECoding/FraminghamAPOE4.dta",  
tab _merge 
tab idtype _merge 
drop _merge 

//Restrict dataset to Offspring Cohort 
keep if idtype==1 
	//n=5,107 Offspring Cohort 
count 
	//n=1,024 w/ MRI 
count if Left_Lateral_Ventricle!=. 
	//n=3,428 w/ APOE
count if APOE!=.
	//n=2,673 w/ DNAm 
count if dunedinpace!=.
count if pchorvath1!=.
	//n=2,379 w/ DNAm & APOE 
count if dunedinpace!=. & APOE!=. 

	//n=948 w/ DNAm & MRI 
count if dunedinpace!=. & Left_Lateral_Ventricle!=. 
	//n=870 w/ DNAm & MRI & APOE
count if dunedinpace!=. & Left_Lateral_Ventricle!=. & APOE!=. 	


//Create variable for age at MRI exam
capture drop mriage 
gen mriage = age1+mri_date/365
label var mriage "Age at MRI"

//Centered age at visit 8 
capture drop cage8
gen cage8 = age8-65

//***************************************************************************//
//Computed MRI variables 
//***************************************************************************//
capture drop MRI 
gen MRI = mri_date!=.

foreach x in TBV HV ICV WMH TSA MCT { 
	capture drop `x'
	}
//Total Brain Volume
gen TBV = TotalGrayVol + CorticalWhiteMatterVol + Left_Cerebellum_White_Matter + Right_Cerebellum_White_Matter
//Intracranial Volume
gen ICV = IntraCranialVol  
//Hippopcampal Volume
gen HV = Left_Hippocampus + Right_Hippocampus
//White Matter Hypointensity 
gen WMH = log(WM_hypointensities)
//Total Surface Area 
#delimit ; 
foreach x in 
	bankssts
	caudalanteriorcingulate
	caudalmiddlefrontal
	cuneus
	entorhinal
	fusiform
	inferiorparietal
	inferiortemporal
	isthmuscingulate
	lateraloccipital
	lateralorbitofrontal
	lingual
	medialorbitofrontal
	middletemporal
	parahippocampal
	paracentral
	parsopercularis
	parsorbitalis
	parstriangularis
	pericalcarine
	postcentral
	posteriorcingulate
	precentral
	precuneus
	rostralanteriorcingulate
	rostralmiddlefrontal
	superiorfrontal
	superiorparietal
	superiortemporal
	supramarginal
	frontalpole
	temporalpole
	transversetemporal
	insula { ; #delimit cr 
	gen wtk_lh_`x' = lh_`x'_tk * lh_`x'_area
	gen wtk_rh_`x' = rh_`x'_tk * rh_`x'_area
	gen area_lh_`x'= lh_`x'_area
	gen area_rh_`x'= rh_`x'_area
	}
egen WTK =rowtotal(wtk_*)
egen TSA = rowtotal(area_*)
gen MCT = WTK / TSA  

label var TBV "Total Brain Volume"
label var HV "Hippocampal Volume"
label var ICV "Intracranial Volume"
label var WMH "White Matter Hypointensities"
label var TSA "Total Surface Area"
label var MCT "Mean Cortical Thickness"

foreach x in TBV HV ICV WMH TSA MCT {
	replace `x'=. if MRI==0
	}

	//Subcortical Volumes 
egen bl_nuc_acc = rowmean(Left_Accumbens_area Right_Accumbens_area)
egen bl_amygdala = rowmean(Left_Amygdala Right_Amygdala)
egen bl = rowmean(Left_Caudate Right_Caudate)
egen bl_cerebellar_cortex = rowmean(Left_Cerebellum_Cortex Right_Cerebellum_Cortex)
egen bl_pallidum = rowmean(Left_Pallidum Right_Pallidum)
egen bl_putamen = rowmean(Left_Putamen Right_Putamen)
egen bl_thalamus = rowmean(Left_Thalamus Right_Thalamus)
egen bl_ventralDC = rowmean(Left_VentralDC Right_VentralDC)
//***************************************************************************//

//CLOCK VARIABLES 
foreach x in pchorvath1 pchannum pcphenoage pcgrimage{
	reg `x' age8
	predict acc_`x', r
	}
gen pace=  dunedinpace 

//Labeled sex variable
capture drop men 
recode sex (2=0), gen(men)
capture label drop men 
label define men 0 "Women" 1 "Men"
label values men men 

//***************************************************************************//
save "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.dta", replace 
//***************************************************************************//


//***************************************************************************//
// EXCEL SHEET RESULTS FOR EW
//***************************************************************************//
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(MRIsample)
putexcel B2 = "The study accession (phs000007.v33.p14) is used to cite the study and its data tables and documents. The data in this file should be cited using the accession pht004364.v2.p14.c1 and pht004364.v2.p14.c2."
putexcel B3 = "Structural Brain Segmentation by MRI Using Freesurfer Software, Offspring Cohort from 2005 - 2010. These include total volume, surface area, gray matter volume, thickness and curvature for the left and right hemispheres. These measures have been assessed by means of two parcellation atlases: Destrieux and Desikan-Killiany."

/*
At MRI exams during 2005-10, we have data for 1,024 individuals
*/
tabstat mriage maxeayrs if idtype ==1 & mri_date!=., s(mean sd n) by(sex) save
putexcel B5 = "MRI Exam 2005-10"
putexcel B6 = "Men"
putexcel B7 = matrix(r(Stat1)), names 
putexcel F6 = "Women"
putexcel F7 = matrix(r(Stat2)), names 
putexcel C7 = "Age at MRI"
putexcel G7 = "Age at MRI"


/*
of the MRI sample, n=948 have DNAm data
*/
tabstat mriage maxeayrs if idtype ==1 & mri_date!=. & pchorvath1!=., s(mean sd n) by(sex) save
putexcel B13 = "MRI Exam w/ DNAm Data"
putexcel B14 = "Men"
putexcel B15 = matrix(r(Stat1)), names 
putexcel F14 = "Women"
putexcel F15 = matrix(r(Stat2)), names 

/*
Compared to the full DNAm sample, the MRI group is younger, more likely to be women, and to have higher education (0.5y), p<0.001 for all
*/


putexcel B21= "MRI+DNAm sample vs DNAm Sample"
putexcel C22 = "B"
putexcel D22 = "SE"
putexcel E22 = "N"
reg age8 MRI if pchorvath1!=. 
	putexcel B23= "Age Diff"
	putexcel C23 = _b[MRI]
	putexcel D23 = _se[MRI]
	local N = e(N)
	putexcel E23 = `N'
reg sex MRI if pchorvath1!=. 
putexcel B24= "Age Diff"
	putexcel C24 = _b[MRI]
	putexcel D24 = _se[MRI]
	local N = e(N)
	putexcel E24 = `N'
reg maxeayrs MRI age8 sex if pchorvath1!=. 
putexcel B25= "Education Diff"
	putexcel C25 = _b[MRI]
	putexcel D25 = _se[MRI]
	local N = e(N)
	putexcel E25 = `N'

	
//Constructed Variables
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(MRIvars)
putexcel B2 = "MRI Variables"
	//All MRI Data
tabstat TBV HV ICV WMH TSA MCT bl_*, s(N mean sd min max) columns(s) save
matrix A = r(StatTotal)'
matrix list A 
putexcel B4=matrix(A), names 
	//DNAm subsample
putexcel B20 = "DNAm + MRI"	
tabstat TBV HV ICV WMH TSA MCT bl_* if pchorvath1!=., s(N mean sd min max) columns(s) save
matrix A = r(StatTotal)'
matrix list A 
putexcel B21=matrix(A), names 

corr TBV HV ICV WMH TSA MCT bl_*
putexcel J4 = matrix(r(C)), names 
corr TBV HV ICV WMH TSA MCT bl_* if pchorvath1!=. 
putexcel J21 = matrix(r(C)), names 
	
//Histograms of MRI Variables 	
foreach y of varlist TBV ICV HV WMH TSA MCT { 
	hist `y', freq scheme(plotplain) name(`y', replace) nodraw 
	}
graph combine TBV ICV HV WMH TSA MCT, scheme(s1mono) cols(2) 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/MRIVarDist.pdf", replace 

//Histograms of Subcortical Volume Variables 	
foreach y of varlist bl_* { 
	hist `y', freq scheme(plotplain) name(`y', replace) nodraw 
	}
graph combine bl_nuc_acc bl_amygdala bl_cerebellar_cortex bl_pallidum bl_putamen bl_thalamus bl_ventralDC, scheme(s1mono) cols(2) 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/MRISubCortVolVarDist.pdf", replace 

//Scatter Plots of MRI Variables by Age 
foreach y of varlist TBV ICV HV WMH TSA MCT { 
	local L: variable label `y'
	#delimit ;
	twoway scatter `y' mriage if sex==2, xtitle(Age at MRI) ytitle(`L') mcolor(pink) msymbol(Oh) 
		|| scatter `y' mriage if sex==1, mcolor(blue) msymbol(+) 
		|| qfit `y' mriage if sex==2, lcolor(cranberry) lpattern(solid) lwidth(medthick) 
		|| qfit `y' mriage if sex==1, lcolor(dknavy) lpattern(solid) lwidth(medthick) scheme(plotplain) legend(off) name(`y', replace) nodraw ; #delimit cr 
	}
graph combine TBV ICV HV WMH TSA MCT, scheme(s1mono) cols(3) xsize(6) ysize(4)
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/MRIVarbyAge.pdf", replace 

//Scatter Plots of Subcortical Volume Variables by Age 
foreach y of varlist bl_* { 
	local L: variable label `y'
	#delimit ;
	twoway scatter `y' mriage if sex==2, xtitle(Age at MRI) ytitle(`L') mcolor(pink) msymbol(Oh) 
		|| scatter `y' mriage if sex==1, mcolor(blue) msymbol(+) 
		|| qfit `y' mriage if sex==2, lcolor(cranberry) lpattern(solid) lwidth(medthick) 
		|| qfit `y' mriage if sex==1, lcolor(dknavy) lpattern(solid) lwidth(medthick) scheme(plotplain) legend(off) name(`y', replace) nodraw ; #delimit cr 
	}
graph combine bl_nuc_acc bl_amygdala bl_cerebellar_cortex bl_pallidum bl_putamen bl_thalamus bl_ventralDC, scheme(s1mono) cols(3) xsize(6) ysize(4)
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/MRISubcortVolVarbyAge.pdf", replace 


//Summary of Clock Data
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(ClockData)

tabstat age8 men pace acc_* if idtype ==1 & pchorvath1!=., s(mean sd n) by(sex) save
putexcel B12 = "DNAm Data"
putexcel B13 = "Men"
putexcel B14 = matrix(r(Stat1)), names 
putexcel F13 = "Women"
putexcel F14 = matrix(r(Stat2)), names 


tabstat age8 men pace acc_* if idtype ==1 & mri_date!=. & pchorvath1!=., s(mean sd n) by(sex) save
putexcel B19 = "DNAm Data + MRI Data"
putexcel B20 = "Men"
putexcel B21 = matrix(r(Stat1)), names 
putexcel F20 = "Women"
putexcel F21 = matrix(r(Stat2)), names 


//Regression analysis
foreach y of varlist TBV HV WMH TSA MCT { 
matrix Fx = J(1,5,999)
	foreach x of varlist pace acc_* {
		capture drop X
		egen X = std(`x') if mri_date!=. & pchorvath1!=. 
		capture drop Y 
		egen Y = std(`y') if mri_date!=. & pchorvath1!=.
		quietly reg Y X c.mriage##c.mriage##men, robust  
		matrix A =  _b[X], _b[X]- invttail(e(df_r),0.025)*_se[X], _b[X]+ invttail(e(df_r),0.025)*_se[X], 2*ttail(e(df_r),abs(_b[X]/_se[X])), e(N)
		matrix rownames A = `x'
		matrix Fx = Fx \ A 
		}
	matrix `y'=Fx[2...,1...]
	matrix colnames `y'=r lb ub p N 
	matrix list `y'
	}
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(MainAnalysis)	
putexcel B1 = "TBV"	
putexcel B2 = matrix(TBV), names	
putexcel B9 = "HV"
putexcel B10 = matrix(HV), names
putexcel B19 = "WMH"
putexcel B20 = matrix(WMH), names
putexcel B29 = "TSA"
putexcel B30 = matrix(TSA), names
putexcel B39 = "MCT"
putexcel B40 = matrix(MCT), names
	
//ADJUSTMENT FOR ICV 
corr TBV ICV
foreach y of varlist TBV HV bl_* { 
matrix Fx = J(1,5,999)
	foreach x of varlist pace acc_* {
		capture drop X
		egen X = std(`x') if mri_date!=. & pchorvath1!=. 
		capture drop Y 
		egen Y = std(`y') if mri_date!=. & pchorvath1!=.
		quietly reg Y X c.mriage##c.mriage##men ICV, robust   
		matrix A =  _b[X], _b[X]- invttail(e(df_r),0.025)*_se[X], _b[X]+ invttail(e(df_r),0.025)*_se[X], 2*ttail(e(df_r),abs(_b[X]/_se[X])), e(N)
		matrix rownames A = `x'
		matrix Fx = Fx \ A 
		}
	matrix `y'=Fx[2...,1...]
	matrix colnames `y'=r lb ub p N 
	matrix list `y'
	}

putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(ICVAdjusted)	
putexcel B1 = "TBV"	
putexcel B2 = matrix(TBV), names	
putexcel B9 = "HV"
putexcel B10 = matrix(HV), names
putexcel B19 = "WMH"
putexcel B20 = matrix(WMH), names
putexcel B29 = "TSA"
putexcel B30 = matrix(TSA), names
putexcel B39 = "MCT"
putexcel B40 = matrix(MCT), names
putexcel B49 = "bl_nuc_acc"	
putexcel B50 = matrix(bl_nuc_acc), names	
putexcel B59 = "bl_amygdala"
putexcel B60 = matrix(bl_amygdala), names
putexcel B69 = "bl_cerebellar_cortex"
putexcel B70 = matrix(bl_cerebellar_cortex), names
putexcel B79 = "bl_pallidum"
putexcel B80 = matrix(bl_pallidum), names
putexcel B89 = "bl_putamen"
putexcel B90 = matrix(bl_putamen), names
putexcel B99 = "bl_thalamus"
putexcel B100 = matrix(bl_thalamus), names
putexcel B109 = "bl_ventralDC"
putexcel B110 = matrix(bl_ventralDC), names

      

//ADJUSTMENT FOR WBC 
foreach y of varlist TBV HV WMH TSA MCT { 
matrix Fx = J(1,5,999)
	foreach x of varlist pace acc_* {
		capture drop X
		egen X = std(`x') if mri_date!=. & pchorvath1!=. 
		capture drop Y 
		egen Y = std(`y') if mri_date!=. & pchorvath1!=.
		quietly reg Y X c.mriage##c.mriage##men  bmem bnv cd4mem cd4nv cd8mem cd8nv nk treg mono bas eos neu 
		matrix A =  _b[X], _b[X]- invttail(e(df_r),0.025)*_se[X], _b[X]+ invttail(e(df_r),0.025)*_se[X], 2*ttail(e(df_r),abs(_b[X]/_se[X])), e(N)
		matrix rownames A = `x'
		matrix Fx = Fx \ A 
		}
	matrix `y'=Fx[2...,1...]
	matrix colnames `y'=r lb ub p N 
	matrix list `y'
	}	
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(WBCAdjusted)	
putexcel B1 = "TBV"	
putexcel B2 = matrix(TBV), names	
putexcel B9 = "HV"
putexcel B10 = matrix(HV), names
putexcel B19 = "WMH"
putexcel B20 = matrix(WMH), names
putexcel B29 = "TSA"
putexcel B30 = matrix(TSA), names
putexcel B39 = "MCT"
putexcel B40 = matrix(MCT), names
	
//ADJUSTMENT FOR APOE4 
tab rs* if MRI==1 & pchorvath1!=. 
tab APOE4 if MRI==1 & pchorvath1!=. 
foreach y of varlist TBV HV WMH TSA MCT { 
matrix Fx = J(1,5,999)
	foreach x of varlist pace acc_* {
		capture drop X
		egen X = std(`x') if mri_date!=. & pchorvath1!=. 
		capture drop Y 
		egen Y = std(`y') if mri_date!=. & pchorvath1!=.
		quietly reg Y X c.mriage##c.mriage##men APOE4 
		matrix A =  _b[X], _b[X]- invttail(e(df_r),0.025)*_se[X], _b[X]+ invttail(e(df_r),0.025)*_se[X], 2*ttail(e(df_r),abs(_b[X]/_se[X])), e(N)
		matrix rownames A = `x'
		matrix Fx = Fx \ A 
		}
	matrix `y'=Fx[2...,1...]
	matrix colnames `y'=r lb ub p N 
	matrix list `y'
	}		
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(APOE4Adjusted)	
putexcel B1 = "TBV"	
putexcel B2 = matrix(TBV), names	
putexcel B9 = "HV"
putexcel B10 = matrix(HV), names
putexcel B19 = "WMH"
putexcel B20 = matrix(WMH), names
putexcel B29 = "TSA"
putexcel B30 = matrix(TSA), names
putexcel B39 = "MCT"
putexcel B40 = matrix(MCT), names
	
//Association with APOE4 
tab rs* if MRI==1 & pchorvath1!=. 
tab APOE4 if MRI==1 & pchorvath1!=. 
matrix Fx = J(1,5,999)
	foreach x of varlist pace acc_* {
		capture drop Y
		egen Y = std(`x') if mri_date!=. & pchorvath1!=.  
		capture drop X 
		gen X = APOE4 
		quietly reg Y X c.mriage##c.mriage##men 
		matrix A =  _b[X], _b[X]- invttail(e(df_r),0.025)*_se[X], _b[X]+ invttail(e(df_r),0.025)*_se[X], 2*ttail(e(df_r),abs(_b[X]/_se[X])), e(N)
		matrix rownames A = `x'
		matrix Fx = Fx \ A 
		}
	matrix Fx=Fx[2...,1...]
	matrix colnames Fx=r lb ub p N 
	matrix list Fx 
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/WhitmanFhamBrainAnalysis.xlsx", modify sheet(APOE4ClockAssoc)	
putexcel B1 = "n=732 w/ 0 E4 alleles, n=140 w/ 1 E4 allele, n=8 w/ 2 E4 alleles"	
putexcel B2 = "rs429358 / rs7412: TT/CC n=622, TT/CT n=99, TT/TT n=1, TC/CC n=130, TC/CT n=10, CC/CC n=8, CC/CT n=0 CC/TT=0"
putexcel B3 = matrix(Fx), names		
	
//Association of APOE4 w/ Dementia in MRI dataset 
stset dem_survdate , failure(dem_status==1) scale(365) origin(time date8)
stcox c.cage8##c.cage8##sex 
stcox APOE4 c.cage8##c.cage8##sex if idtype==1 & pchorvath1!=. & MRI==1

	
//Visit Dates of MRI and DNAm collection 
scatter mri_date date8, msymbol(O) mcolor(navy) || lfit date8 date8, legend(off) lwidth(thick) lpattern(solid) lcolor(cranberry) xtitle("Date of Visit 8" "(days since 1st visit)")  ytitle("Date of MRI" "(days since 1st visit)") 	scheme(plotplain)  
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/EthanWhitman/FraminghamBrainImaging/R2/MRIDatebyVisit8Date.pdf", replace 
	


