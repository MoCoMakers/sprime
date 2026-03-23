import pandas as pd

path = "excel ref.xlsx"
sheet = "Sprime_delta_Sprime"

df = pd.read_excel(path, sheet_name=sheet, header=None)

keywords = ["IC50", "4-PL", "Parameter", "Values", "S' for", "delta", "ΔS"]

hits = []
for i in range(len(df)):
    row = df.iloc[i]

    
    row_text = " | ".join(row.dropna().astype(str)).lower()

    if any(k.lower() in row_text for k in keywords):
        hits.append(i)

print("Candidate header/info rows:", hits[:50])

for i in hits[:10]:
    print("\n--- around row", i, "---")
    print(df.iloc[max(0, i-2): i+3, :25])


HEADER_ROW = 5   # ← 여기 바꿔보면서 맞추기

df = pd.read_excel(path, header=HEADER_ROW)
print(df.head())
