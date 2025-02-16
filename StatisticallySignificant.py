import pandas as pd
import numpy as np
import scipy.stats as stats
# Function to calculate Odds Ratio and 95% Confidence Interval
def odds_ratio_and_ci(a, b, c, d, confidence_level=0.95):
   # Continuity correction if any value is zero
   if a == 0 or b == 0 or c == 0 or d == 0:
       a += 0.5
       b += 0.5
       c += 0.5
       d += 0.5
   # Odds Ratio
   odds_ratio = (a * d) / (b * c) if (b * c) != 0 else np.nan
   # Standard error of the log odds ratio
   se = np.sqrt((1/a) + (1/b) + (1/c) + (1/d))
   # Z-score for confidence interval
   z = stats.norm.ppf(1 - (1 - confidence_level) / 2)
   # Confidence Interval
   ci_lower = np.exp(np.log(odds_ratio) - z * se) if odds_ratio > 0 else np.nan
   ci_upper = np.exp(np.log(odds_ratio) + z * se) if odds_ratio > 0 else np.nan
   return odds_ratio, (ci_lower, ci_upper)
# Load dataset (expects an Excel file with columns: Gene, North_Present, North_Absent, South_Present, South_Absent)
file_path = 'contingencytablechapter3.xlsx'
sheet_name = 'Sheet1'
df = pd.read_excel(file_path, sheet_name=sheet_name)
def analyze_genes(df):
   results = []
   for index, row in df.iterrows():
       print(f"Processing Gene: {row['Gene']}")
       # Contingency Table for North Island vs South Island
       table = np.array([
           [row['North_Present'], row['North_Absent']],
           [row['South_Present'], row['South_Absent']]
       ])
       try:
           # Perform Chi-Square Test
           chi2_stat, p_value_chi2, _, expected = stats.chi2_contingency(table, correction=False)
           # If expected values are < 5, use Fisher’s Exact Test
           if np.any(expected < 5):
               test_name = "Fisher's Exact"
               oddsratio_fisher, p_value_fisher = stats.fisher_exact(table, alternative='two-sided')
               p_value = p_value_fisher
               chi2_stat = np.nan  # Chi-square isn't valid for small expected counts
           else:
               test_name = "Chi-Square"
               p_value = p_value_chi2
               p_value_fisher = np.nan  # Fisher's test not needed if Chi-square is valid
       except Exception as e:
           print(f"Exception encountered: {e}")
           test_name = "Fisher's Exact"
           oddsratio_fisher, p_value_fisher = stats.fisher_exact(table, alternative='two-sided')
           p_value = p_value_fisher
           chi2_stat = np.nan
           p_value_chi2 = np.nan
       # Calculate Odds Ratio and 95% Confidence Interval
       a, b, c, d = table.flatten()
       odds_ratio, (ci_lower, ci_upper) = odds_ratio_and_ci(a, b, c, d)
       # Determine significance
       is_significant = p_value < 0.05
       # Store results
       results.append({
           'Gene': row['Gene'],
           'Test Used': test_name,
           'P-Value': p_value,
           'Chi-Square Statistic': chi2_stat,
           'Chi-Square P-Value': p_value_chi2,
           "Fisher's Exact P-Value": p_value_fisher,
           'Odds Ratio': odds_ratio,
           '95% CI Lower': ci_lower,
           '95% CI Upper': ci_upper,
           'Significance': is_significant,
           'Expected Counts': expected.tolist()
       })
   # Convert results to DataFrame
   results_df = pd.DataFrame(results)
   return results_df
# Analyze data and save results
results_df = analyze_genes(df)
print("Results DataFrame:")
print(results_df)
# Save results to CSV
results_df.to_csv('Gene_Comparison_Results.csv', index=False)
