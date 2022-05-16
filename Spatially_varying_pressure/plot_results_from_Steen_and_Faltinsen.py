import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('C:/Users/mathi/OneDrive - NTNU/Master Thesis/Spatially varying pressure/Results from Steen and Faltinsen (1995)/Rigid panel model/Uniform pressure/uniform_pressure_RAO_nice_format.csv')

print(df.iloc[:, 1])

plt.plot(df.iloc[:, 1], df.iloc[:, 2])
plt.title('RAO from Steen and Faltinsen (1995)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$\\mu_{ua}$ [-] / $\\zeta_a$ [m]')
plt.show()



