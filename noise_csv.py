import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('output.csv', header=None)
df.columns = ['hpz', 'hxz', 'hpx', 'hxx']

dt = 0.0005
time = [i * dt for i in range(int(len(df)))]

data1 = df['hpz']
data2 = df['hxz']
data3 = df['hpx']
data4 = df['hxx']

#print(data1, data2, data3, data4)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.plot(time, data3, data1, label='h_plus_combined')
ax.plot(time, data4, data2, label='h_x_combined')
ax.plot(time, 0, data1, alpha = 0.4, label='h_plus_z-axis')
ax.plot(time, 0, data2, alpha = 0.4, label='h_x_z-axis')
ax.plot(time, data3, 0, alpha = 0.4, label='h_plus_x-axis')
ax.plot(time, data4, 0, alpha = 0.4, label='h_x_x-axis')

ax.set_xlabel('t[years]')
ax.set_ylabel('X')
ax.set_zlabel('Z')

ax.view_init(elev=20., azim=-35, roll=0)
ax.legend()

plt.grid(True)
plt.show()