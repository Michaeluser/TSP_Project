import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data_frame = pd.read_csv("SA_succession.txt")

generation_column = list(data_frame["Generation"]) #x variable
length_column = list(data_frame["Length"]) #y variable
fitness_column = np.array([int(10**12 / length_column[gen]) for gen in range(len(generation_column))])

Q1 = np.percentile(fitness_column, 25)
Q3 = np.percentile(fitness_column, 75)
IQR = Q3 - Q1

l_bound = Q1 - 1.5 * IQR
u_bound = Q3 + 1.5 * IQR

outliers = fitness_column[(fitness_column > u_bound) | (fitness_column < l_bound)]

print(len(outliers))

plt.plot(generation_column, fitness_column)
plt.ylabel("Generation fitness")
plt.xlabel("Generation number")
plt.show()