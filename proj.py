#immune infiltration signatures predict survival in breat cancer - we will look at the expression of immune genes in mrna data, and survival times (DSS_MONTHS and OS_STATUS) in clinical data
#immune infiltration just means expression of immune genes


#do patients with lots of immune cells/expression live longer?
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

#==================== LOAD DATA (GATHERED FROM CBIOPORTAL) ==========================

#expression data
exp = pd.read_csv(
    r"C:/Users/xavr/Documents/brca/data/data_mrna_seq_v2_rsem.txt",
    sep="\t",
    comment="#"
)

#clinical data 
clin = pd.read_csv(
    r"C:/Users/xavr/Documents/brca/data/data_clinical_patient.txt",
    sep="\t",
    comment="#"
)


#for mac

#exp = pd.read_csv(r"/Users/user/Desktop/Bioinformatics/brca/data_mrna_seq_v2_rsem.txt",
#    sep="\t",
#    comment="#"
#)
#clin = pd.read_csv(r"/Users/user/Desktop/Bioinformatics/brca/data_clinicaL_patient.txt",
#    sep="\t",
#    comment="#"
#)

#==================== EXTRACT SURVIVAL DATA ==========================

#we will use 3 survival variables from the patient data for survival analysis
    #1. DSS_MONTHS, records time from initial diagnosis to time of cancer-specific death
    #2. OS_STATUS, records patient death due to cancer and potentially non-cancer specific causes
    #3. DSS_STATUS, records cancer-related death, but is ambiguous as it gives a dead or alive status for many patients. This will be used together with OS_STATUS to confirm cancer-related death 

#first check for missing values from our survival data (DSS_MONTHS, OS_STATUS, DSS_STATUS)
print("CHECK FOR MISSING VALUES")

print(f"DSS_MONTHS missing: {clin['DSS_MONTHS'].isna().sum()}")
print(f"OS_STATUS missing: {clin['OS_STATUS'].isna().sum()}")
print(f"DSS_STATUS missing: {clin['DSS_STATUS'].isna().sum()}") #MISSING 20 PATIENTS

#using .dropna, remove the 20 out of 1084 patients missing DSS_STATUS values
clin_cleaned = clin.dropna(subset = ["DSS_STATUS"]).copy()


#months from initial diagnosis to cancer-related death
months = clin_cleaned["DSS_MONTHS"]


#composite event variable function, where event is cancer death
def cancer_death(row):
    os_s = int(str(row["OS_STATUS"]).split(":")[0])
    dss = int(str(row["OS_STATUS"]).split(":")[0])

    #if os_s = 1("DECEASED") and dss = 1 ("DEAD WITH TUMOR")
    if os_s == 1 and dss == 1:
        return 1
    else:
        return 0


#create new column "tumor_death" in cleaned dataset using the function
clin_cleaned["tumor_death"] = clin_cleaned.apply(cancer_death, axis=1)
#print(clin_cleaned["tumor_death"].value_counts(dropna=False))












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