#! /usr/bin/env python2.7

import matplotlib.pyplot as plt
import pandas as pd

df  = pd.read_csv("bin/solution", names=["solution_value"])
idx = [float(x)/len(df) for x in range(1,len(df)+1)]

df.index = idx
df.index.name = "x"

df.plot(style=".-", figsize=(5,3))
plt.savefig("plot/solution.png")

