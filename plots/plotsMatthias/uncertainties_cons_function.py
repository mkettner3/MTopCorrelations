import matplotlib.pyplot as plt
import numpy as np

plt.figure()
for v, c in zip([-2, -1, -0.5, 0.5, 1, 2], ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']):
    plt.plot(np.linspace(0, 10, num=100), [1+np.sign(v) * 0.01]*100, color=c, label=f'v = {v}')
    plt.plot(np.linspace(10, 30, num=200), 1+np.sign(v) * (0.01 + 0.03/90 * (np.linspace(10, 30, num=200)-10) * np.abs(v)), color=c)
plt.xlabel('$p_T$')
plt.ylabel('k')
plt.legend()
plt.show()
