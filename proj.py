#immune infiltration signatures predict survival in breat cancer - we will look at the expression of immune genes in mrna data, and survival times (OS_MONTHS and OS_STATUS) in clinical data
#immune infiltration just means expression of immune genes


#do patients with lots of immune cells/expression live longer?
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

##load data from cBioPortal

#expression data
#exp = pd.read_csv(
#    r"C:/Users/xavr/Documents/brca/data/data_mrna_seq_v2_rsem.txt",
#    sep="\t",
#    comment="#"
#)

#clinical data 
#clin = pd.read_csv(
#    r"C:/Users/xavr/Documents/brca/data/data_clinical_patient.txt",
#    sep="\t",
#    comment="#"
#)


#for mac

exp = pd.read_csv(r"/Users/user/Desktop/Bioinformatics/brca/data_mrna_seq_v2_rsem.txt",
    sep="\t",
    comment="#"
)
clin = pd.read_csv(r"/Users/user/Desktop/Bioinformatics/brca/data_clinicaL_patient.txt",
    sep="\t",
    comment="#"
)
print(clin["OS_MONTHS"])

#extract survival data

#first check for missing values from our survival data (OS_MONTHS, OS_STATUS, DSS_STATUS)

#check 
print("CHECK FOR MISSING VALUES")

print(f"OS_MONTHS missing: {clin['OS_MONTHS'].isna().sum()}")
print(f"OS_STATUS missing: {clin['OS_STATUS'].isna().sum()}")


print(f"DSS_STATUS missing: {clin['DSS_STATUS'].isna().sum()}")

#using .dropna, remove the 20 out of 1084 patients missing DSS_STATUS values
clin_cleaned = clin.dropna(subset = ["DSS_STATUS"])


#months from initial diagnosis to recorded death
months = clin_cleaned["OS_MONTHS"]

#extract dead or alive survival status (includes non-cancer related deaths)
os_stat = clin_cleaned["OS_STATUS"] #0 = ALIVE, 1 = DEAD
#simplify OS_STAT status to just the 0s and 1s
os_stat_val = clin_cleaned["OS_STATUS"].str.split(':').str[0].astype(int)

#patient's dead or alive status (includes only cancer-related deaths, but needs to be tied to os_stat due for ambiguous "dead/alive tumor free" patient status)
dss_stat = clin_cleaned["DSS_STATUS"]
#simplify DSS_STATUS to 0s and 1s
dss_stat_val = dss_stat.str.split(':').str[0].astype(int)



#composite event variable function
#this is the variable we will use to confirm death by cancer

def cancer_death(row):


    if os_stat_val and dss_stat_val == 0:
        return 0 #ALIVE, cancer free
    elif os_stat_val == 0 and dss_stat_val == 1:
        return 0 #alive, cancer free
    elif os_stat_val and dss_stat_val == 1:
        return 1 #dead, had cancer
    #check data to make sure no rows where OS_STATUS= 0 and DSS_STATS= 1 (contradiction, alive and dead with tumor)

#os_stat_val
#0:LIVINg 1:DEAD WITH TUMOR

#dsstat_val0:ALIVE OR DEAD TUMOR FREE1:DECEASED








# notebook 2: Calculate immune scores
# Use ESTIMATE algorithm or immunedeconv package
# This gives you "immune score" for each patient

# notebook 3: Survival analysis
# Split patients into High vs Low immune groups
# Plot Kaplan-Meier curves
# Run Cox regression to check significance

# notebook 4: Make beautiful figures
# Heatmap of immune genes
# Boxplots by cancer subtype
# Final summary figure

print("="*50)