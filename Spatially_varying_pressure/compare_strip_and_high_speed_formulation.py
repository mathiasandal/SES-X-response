import pandas as pd
import matplotlib.pyplot as plt

df_strip_theory = pd.read_csv('Results/results_strip_theory.csv')
df_high_speed = pd.read_csv('Results/results_high_speed.csv')



plt.plot(df_strip_theory.iloc[:, 1], df_strip_theory.iloc[:, 3], label='strip')
plt.plot(df_high_speed.iloc[:, 1], df_high_speed.iloc[:, 3], label='high speed')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('eta_3a [m]')
plt.show()


plt.plot(df_strip_theory.iloc[:, 1], df_strip_theory.iloc[:, 4], label='strip')
plt.plot(df_high_speed.iloc[:, 1], df_high_speed.iloc[:, 4], label='high speed')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('eta_5a [rad]')
plt.show()

plt.plot(df_strip_theory.iloc[:, 1], df_strip_theory.iloc[:, 5], label='strip')
plt.plot(df_high_speed.iloc[:, 1], df_high_speed.iloc[:, 5], label='high speed')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('mu_ua [-]')
plt.show()