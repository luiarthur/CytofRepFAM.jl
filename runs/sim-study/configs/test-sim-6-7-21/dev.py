import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

df_stacked = pd.read_csv("img/df.csv")

sns.boxplot(data=df_stacked, x="feature", y="y", hue="marker")
plt.axhline(0)
plt.savefig("img/express.pdf", bbox_inches="tight")
plt.close()
