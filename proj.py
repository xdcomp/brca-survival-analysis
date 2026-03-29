#immune infiltration signatures predict survival in breat cancer - we will look at the expression of immune genes in mrna data, and survival times (DSS_MONTHS and OS_STATUS) in clinical data
#immune infiltration just means expression of immune genes


#do patients with lots of immune cells/expression live longer?
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

#==================== LOAD DATA (GATHERED FROM CBIOPORTAL) ==========================

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


#for mac, make sure to move to 

exp = pd.read_csv(r"/Users/user/Desktop/Bioinformatics/brca/data_mrna_seq_v2_rsem.txt",
    sep="\t",
    comment="#"
)
clin = pd.read_csv(r"/Users/user/Desktop/Bioinformatics/brca/data_clinicaL_patient.txt",
    sep="\t",
    comment="#"
)

#==================== EXTRACT SURVIVAL DATA ==========================

#we will use 3 survival variables from the patient data for survival analysis
    #1. DSS_MONTHS - records time from initial diagnosis to time of cancer-specific death
    #2. OS_STATUS - clear dead or alive status, with death to cancer and potentially non-cancer specific causes
    #3. DSS_STATUS - either clear death with tumor or dead/alive status with tumor. This will be used together with OS_STATUS to confirm cancer-related death 

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


#==================== CALCULATE IMMUNE SCORES ==========================


#we are using the exp data (patients are columns, genes are rows. we will transpose in later step, as analysis tools expect patients as rows)
##print(exp.head())




#create list of immune genes to look for in the data
immune_genes = ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG",
                  "FOXP3", "IL2RA", "CTLA4",
                  "CD19", "MS4A1", "CD79A",
                  "NCAM1", "NKG7", "KLRD1",
                  "CD68", "CD163", "MRC1"]

#check what genes are avaiable in the data
available = exp["Hugo_Symbol"].isin(immune_genes) #all genes are available!

#print("here's what's AVAILABLE:\n\n\n")
print(exp[available]["Hugo_Symbol"].tolist())

#cell types dictionary
immune_signatures = {
    "CD8_T_cells":  ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG"],
    "Tregs":        ["FOXP3", "IL2RA", "CTLA4"],
    "B_cells":      ["CD19", "MS4A1", "CD79A"],
    "NK_cells":     ["NCAM1", "NKG7", "KLRD1"],
    "Macrophages":  ["CD68", "CD163", "MRC1"]
}

#filters exp dataframe down to only the 17 immune genes
exp_immune = exp[exp["Hugo_Symbol"].isin(immune_genes)]


#drop the Entrez_Gene_Id column, set gene names/Hugo_Symbol to row labels # now there are only gene names, patient ids, and genes
exp_immune = exp_immune.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"])


#transpose the data so analysis tools may use the patient rows 
exp_immune = exp_immune.T

#print(exp_immune.head())
#print(exp_immune.shape) # gives (1082,17). the original exp data i downloaded only contained 1082, might have to lookup how their data was lost in the cBioPortal webpage for data


#log-transform the data ()
import numpy as np

# the +1 is to avoid log2(0) which gives negative infinity
exp_immune_log = np.log2(exp_immune + 1)

#calculate average immune scores per cell type

#for each cell type in cell types dictionary...
for cell_type, genes in immune_signatures.items():# items = genes

    #create new column for each cell type, with row value being the averaged gene expression of iterated cell type
    exp_immune_log[cell_type] = exp_immune_log[genes].mean(axis=1) 

#pull keys (cell types) and convert to a list
score_cols = list(immune_signatures.keys())

#final dataframe, subset into the 5 existing cell_type score columns
immune_scores = exp_immune_log[score_cols]

#preview first 5 rows of scores
print("printing immune SCORES")
print(immune_scores.head())
print("DONE PRINTING IMMUNE SCORES")

#========================= MERGE THE DATA ===============================

#before merging, we confirm patient IDs in clin_cleaned match ID format in immune_scores
immune_scores.index = immune_scores.index.str.replace('-01', '', regex=False)

#now we can merge
merged = clin_cleaned.merge(immune_scores, left_on="PATIENT_ID", right_index=True, how="inner")
print("printing merged shape")
print(merged.shape)
print("done")

print("printing merged head")
print(merged.head())
print("done")


print("PRINTING MERGE DATA")
print(merged[["PATIENT_ID", "OS_MONTHS", "tumor_death", "CD8_T_cells", "Tregs", "B_cells", "NK_cells", "Macrophages"]].head())
print("done thanks")



#========================= SURVIVAL ANALYSIS ===============================

# notebook 3: Survival analysis
# Split patients into High vs Low immune groups FOR EACH CELL TYPE
# Plot Kaplan-Meier curves
# Run Cox regression to check significance

# notebook 4: Make beautiful figures
# Heatmap of immune genes
# Boxplots by cancer subtype
# Final summary figure

print("="*50)

print(clin_cleaned["PATIENT_ID"].head())
print(immune_scores.index[:5]) #looks like there is a mismatch in expression and clinical formatting of patient data "Clinical data has TCGA-3C-AAAU but expression data has TCGA-3C-AAAU-01."