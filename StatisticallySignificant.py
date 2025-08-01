import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

# Function to calculate odds ratio and confidence interval
def odds_ratio_and_ci(a, b, c, d, confidence_level=0.95):
    # Apply continuity correction if any count is zero
    if a == 0 or b == 0 or c == 0 or d == 0:
        a += 0.5
        b += 0.5
        c += 0.5
        d += 0.5
    odds_ratio = (a * d) / (b * c)
    se = np.sqrt((1/a) + (1/b) + (1/c) + (1/d))
    z = stats.norm.ppf(1 - (1 - confidence_level) / 2)
    ci_lower = np.exp(np.log(odds_ratio) - z * se)
    ci_upper = np.exp(np.log(odds_ratio) + z * se)
    return odds_ratio, ci_lower, ci_upper

# Load data from Excel
file_path = 'contingencytablechapter3.xlsx'
sheet_name = 'Sheet2'
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Analyze genes
results = []
raw_p_values = []

for index, row in df.iterrows():
    gene = row['Gene']
    a = row['North_Present']
    b = row['North_Absent']
    c = row['South_Present']
    d = row['South_Absent']
    
    table = np.array([[a, b], [c, d]])
    
    try:
        chi2_stat, chi2_p, _, expected = stats.chi2_contingency(table, correction=False)
        if np.any(expected < 5):
            test_used = "Fisher's Exact"
            oddsratio, fisher_p = stats.fisher_exact(table)
            p_value = fisher_p
            chi2_stat = np.nan
            chi2_p = np.nan
        else:
            test_used = "Chi-Square"
            oddsratio, fisher_p = np.nan, np.nan
            p_value = chi2_p
    except Exception as e:
        print(f"Error processing gene {gene}: {e}")
        test_used = "Fisher's Exact"
        oddsratio, fisher_p = stats.fisher_exact(table)
        p_value = fisher_p
        chi2_stat = np.nan
        chi2_p = np.nan
        expected = np.full((2, 2), np.nan)
    
    or_val, ci_lower, ci_upper = odds_ratio_and_ci(a, b, c, d)
    raw_p_values.append(p_value)
    
    results.append({
        'Gene': gene,
        'Test Used': test_used,
        'P-Value': p_value,
        'Chi-Square Statistic': chi2_stat,
        'Chi-Square P-Value': chi2_p,
        "Fisher's Exact P-Value": fisher_p,
        'Odds Ratio': or_val,
        '95% CI Lower': ci_lower,
        '95% CI Upper': ci_upper,
        'Expected Counts': expected.tolist() if isinstance(expected, np.ndarray) else None
    })

# FDR Correction
pvals = np.array(raw_p_values)
_, fdr_pvals, _, _ = multipletests(pvals, method='fdr_bh')

# Add FDR-adjusted p-values and significance
for i, result in enumerate(results):
    result['FDR Adjusted P-Value'] = fdr_pvals[i]
    result['Significance (FDR < 0.05)'] = fdr_pvals[i] < 0.05

# Convert to DataFrame and save
results_df = pd.DataFrame(results)
print(results_df)
results_df.to_csv('Gene_Comparison_Results_teat_spray.csv', index=False)
