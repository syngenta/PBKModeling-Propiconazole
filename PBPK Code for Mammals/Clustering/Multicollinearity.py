import sklearn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import os
import math
from functools import partial, reduce
from imblearn.over_sampling import SMOTE
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, r2_score
import statsmodels
from statsmodels.stats.outliers_influence import variance_inflation_factor

#####################    Data preprocessing   ########################

##############################################
##              Data loading                ##
##############################################
# dataset path
cwd = 'C:\\Users\\XXX\OneDrive - Syngenta\\AI\\Propiconazole\\PBK_Modeling\\in vitro'
dat_path = cwd + '\\4clustering_1007.csv'

# read data
dat = pd.read_csv(dat_path, encoding = 'unicode_escape', engine='python')
X   = dat.drop(columns=['Compound'])    # Conc_xxx are target variables

#X   = dat.drop(columns=['Compound', 'Acute_LD50', 'Chronic_NOAEL'])
ax=plt.figure(figsize=(12,15))
corr_matrix = X.corr()
# Print the correlation matrix
print(corr_matrix)

# Optionally, you can visualize the correlation matrix using a heatmap
rdgn = sns.diverging_palette(h_neg=240, h_pos=100, s=99, l=55, sep=3, as_cmap=True)
sns.heatmap(corr_matrix , linewidth=0.5, cmap=rdgn, center=0.00, annot=True, fmt ='.0%', annot_kws={"size": 8})
plt.subplots_adjust(left=0.2, bottom=0.25)
plt.xticks(rotation = 90, fontsize=12)
plt.yticks(rotation = 0, fontsize=12)
plt.show()
ax.savefig('C:\\Users\\XXX\OneDrive - Syngenta\\AI\\Propiconazole\\PBK_Modeling\\Read-cross code\\Heatmap_overall.tiff', dpi = 600, pil_kwargs={"compression": "tiff_lzw"})


df = X

def drop_highly_correlated(df, threshold):
    # Create correlation matrix
    corr_matrix = df.corr().abs()

    # Select upper triangle of correlation matrix
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool_))

    # Find index of feature columns with correlation greater than threshold
    to_drop = [column for column in upper.columns if any(upper[column] > threshold)]

    # Drop features 
    df.drop(to_drop, axis=1, inplace=True)

    return df

# Usage:
X= drop_highly_correlated(X, 0.95)

ax=plt.figure(figsize=(12,15))
corr_matrix = X.corr()
# Print the correlation matrix
print(corr_matrix)

# Optionally, you can visualize the correlation matrix using a heatmap
rdgn = sns.diverging_palette(h_neg=240, h_pos=100, s=99, l=55, sep=3, as_cmap=True)
sns.heatmap(corr_matrix , linewidth=0.5, cmap=rdgn, center=0.00, annot=True, fmt ='.0%')
plt.subplots_adjust(left=0.2, bottom=0.25)
plt.xticks(rotation = 90, fontsize=12)
plt.yticks(rotation = 0, fontsize=12)
plt.show()
ax.savefig('C:\\Users\\XXX\OneDrive - Syngenta\\AI\\Propiconazole\\PBK_Modeling\\Read-cross code\\Heatmap_filtering.tiff', dpi = 600, pil_kwargs={"compression": "tiff_lzw"})


df = X
# Calculate correlation matrix
corr_matrix = df.corr()
# Find pairs of perfectly correlated columns
perfect_corr_pairs = np.where(np.isclose(corr_matrix, 1.0))
# Get the variable names
variables = [(corr_matrix.columns[x], corr_matrix.columns[y]) for x, y in zip(*perfect_corr_pairs) if x != y and x < y]
# Print the variable names
for var1, var2 in variables:
    print(f"{var1} and {var2} are perfectly correlated.")

'''BDDCS  and ECCS are perfectly correlated.
BDDCS  and observed_clearance_mechanism are perfectly correlated.
ECCS and observed_clearance_mechanism are perfectly correlated.'''
    

# check if any column has types as object
object_columns = X.select_dtypes(include=['object']).columns
print(object_columns)

'''MW	Solubility (mol/L)	LogP	pKa_Accept	Vapor pressure	fup_human	Permeability	Clint_human	Taminoto_Score 	Acute_LD50 	Chronic_NOAEL
'''

###################################################
##              Multicollinearity                ##
###################################################
# reference: https://github.com/siddharthfcb/WSOM_RF/blob/main/Random_forest_GPR_SVR.ipynb
# https://towardsdatascience.com/hyperparameter-tuning-the-random-forest-in-python-using-scikit-learn-28d2aa77dd74
# https://www.datacareer.de/blog/random-forest-in-python-with-scikit-learn/

def calc_vif(X):

    # Calculating VIF
    vif = pd.DataFrame()
    vif["variables"] = X.columns
    vif["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]

    return(vif)

vif = calc_vif(X)
print(vif)
vif.to_csv('C:\\Users\\XXX\OneDrive - Syngenta\\AI\\Propiconazole\\PBK_Modeling\\Read-cross code\\vif_notox.csv', index=False)

# Plotting VIF values
plt.figure(figsize=(10, 6))
plt.bar(vif['variables'], vif['VIF'])
plt.xlabel('Variables')
plt.ylabel('VIF')
plt.title('VIF values')
plt.xticks(rotation=90)
plt.show()




# drop varibale with high VIF
dropped_variables = pd.DataFrame(columns=['variable', 'VIF'])

while True:
    vif     = calc_vif(X)
    max_vif = vif['VIF'].max()
    if max_vif > 10:
        max_vif_var       = vif[vif['VIF'] == max_vif]['variables'].values[0]
        new_row           = {'variable': max_vif_var, 'VIF': max_vif}
        dropped_variables = pd.concat([dropped_variables, pd.DataFrame([new_row])], ignore_index=True, axis=0)
        X                 = X.drop(columns=max_vif_var)
        print("Dropped variables:", max_vif_var)
    else:
        print("=======   END=========")
        break

print("Dropped variables:")
print(dropped_variables)
print("\nFinal VIFs:")
print(vif)
dropped_variables.to_csv('C:\\Users\\XXX\OneDrive - Syngenta\\AI\\Propiconazole\\PBK_Modeling\\Read-cross code\\dropped_variables_wTOX.csv', index=False)
vif.to_csv('C:\\Users\\XXX\OneDrive - Syngenta\\AI\\Propiconazole\\PBK_Modeling\\Read-cross code\\vif_final_wTOX.csv', index=False)
