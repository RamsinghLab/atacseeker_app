import pandas as pd 
import sys

if len(sys.argv) < 2:
	print "Usage: python" + sys.argv[0] + "S_file N_file"


df1 = pd.read_table(sys.argv[1], sep = ",")
df2.columns = ["S-Species", "S-Individual", "Position", "S-Major", "S-Minor", 
	"S-Ratio", "S-Coverage", "S-Bases"]

df2 = pd.read_table(sys.argv[2], sep = ",")
df2.columns = ["N-Species", "N-Individual", "Position", "N-Major", "N-Minor", 
	"N-Ratio", "N-Coverage", "N-Bases"]

df = pd.merge(df1, df2, on = "Position", how = "outer")
df.fillna("-", inplace = True)

df[["Position", "S-Major", "S-Minor", "S-Ratio", 
	"N-Major", "N-Minor", "N-Ratio", 
	"S-Coverage", "S-Bases", 
	"N-Coverage", "N-Bases"]].to_csv("het.out", index = False, sep "\t")
