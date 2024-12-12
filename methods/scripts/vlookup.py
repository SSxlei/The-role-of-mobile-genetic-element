import pandas as pd

file_path = 'Seep.xlsx'
try:
    df = pd.read_csv(file_path, delimiter='\t', error_bad_lines=False, warn_bad_lines=True)
except pd.errors.ParserError as e:
    print(f"Error parsing the CSV file: {e}")

with open('MGEid', 'r') as file:
    mgeid_list = [line.strip() for line in file]

results = []

for mgeid in mgeid_list:
    try:
        matched_rows = df[df.iloc[:, 0] == mgeid]
        if not matched_rows.empty:
            for _, row in matched_rows.iterrows():
                results.append(f"{mgeid}\t{row.iloc[6]}")
        else:
            results.append(f"{mgeid}\tNone")
    except Exception as e:
        print(f"Error processing MGEid {mgeid}: {e}")

with open('MGE.txt', 'w') as output_file:
    for result in results:
        output_file.write(result + '\n')
