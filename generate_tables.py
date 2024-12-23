import pandas as pd
import numpy as np

def normalize_name(name: str):
    return name.split("/")[1].split(".")[0]

df = pd.DataFrame()
times_df = pd.DataFrame()
ks = [2**i for i in range(9, 10)]

for k in ks:
    # Load the CSV files
    results_file = f'results/synthetic-{k}.csv'  # Replace with your results file path
    times_file = f'results/synthetic-times-{k}.csv'  # Replace with your times file path

    # Read the data into Pandas dataframes
    results_df = pd.read_csv(results_file)
    #results_df["k"] = k
    df = pd.concat([df, results_df])
    aux_times_df = pd.read_csv(times_file)
    #aux_times_df["k"] = k
    times_df = pd.concat([times_df, aux_times_df])

#df['input_file'] = df["input_file"].apply(normalize_name)
#times_df['input_file'] = times_df["input_file"].apply(normalize_name)

# Group data by `input_file` and calculate mean and SEM for results and times
cols = [c for c in df.columns if c not in ["n", " N"]]
print(cols)
#for col in cols:
 #   df[col] = ((df[col] - df[" Card"])**2)
results_stats = df.groupby(['n', " N"]).mean().reset_index()
"""
results_stats[cols] = ((results_stats[cols]/100)**0.5) # SE
for col in cols:
    results_stats[col] = results_stats[col]/results_stats[" Card"]
results_stats = results_stats.drop(columns=[" Card"])
"""
results_stats = results_stats.groupby(["n", " N"]).mean()

times_stats = times_df.groupby(['n', " N"]).mean()

# Generate LaTeX table from the combined statistics
latex_table = results_stats.to_latex(float_format="%.0f", index_names=False)

# Save the LaTeX table to a file
with open("statistics_table.tex", "w") as f:
    f.write(latex_table)

print("LaTeX table saved to 'statistics_table.tex'")
