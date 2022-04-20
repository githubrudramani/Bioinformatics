import pandas as pd
import glob
print("Merging started")
def merge_count(input_files = "*txt"):
    files = glob.glob(input_files)
    L = []
    for file in files:
        print(f"reading file: {file}")
        try:
            f = pd.read_csv(file,sep = "\s+", header = None, index_col = 0)
            f.columns = [file]
            L.append(f)
        except:
            pass
    df = pd.concat(L)
    df = df.fillna(0)
    discard = []
    for i in df.index:
        if  i.startswith("_"):
            discard.append(i)
    df = df.drop(discard)
    df =df.groupby(df.index).mean()
    df.columns = [i[0:-4] for i in df.columns]
    return df
df = merge_count("*txt")
df.to_csv("merged_counts.csv")
meta = pd.DataFrame(df.columns, columns = ["samples"])
meta.to_csv("metaData.csv")
print("Data merged succesfully")
