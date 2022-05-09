import pandas as pd
import matplotlib.pyplot as plt

filepath = 'C:/Users/mathi/OneDrive - NTNU/Master Thesis/Spatially varying pressure/Results from Steen and Faltinsen ' \
           '(1995)/Rigid panel model/Uniform pressure/'
filename = 'uniform_pressure_RAO.csv'

df = pd.read_csv(filepath + filename)

print(df)

new_df = pd.DataFrame()
new_df['Encounter frequency [Hz]'] = df['x'].apply(lambda x: float(x.replace(',', '.')))
new_df['mu_u [-] / wave amp. [m]'] = df['Uniform_pressure_RAO'].apply(lambda x: float(x.replace(',', '.')))

new_df.to_csv(filepath + 'uniform_pressure_RAO_nice_format.csv')

plt.plot(new_df.iloc[:, 0], new_df.iloc[:, 1], 'kx-')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$mu_a$ [-] / $\\zeta_a$ [m]')
plt.xlim([new_df.iloc[0, 0], new_df.iloc[-1, 0]])
plt.show()
