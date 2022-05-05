# -*- coding: utf-8 -*-

# Python 3 script

import numpy as np
import matplotlib.pyplot as plt

plt.figure()
plt.plot([1.984, 1.804, 1.673, 1.559, 1.455, 1.404])
plt.ylabel('Mean Peak Position of $3\zeta$')
plt.xlabel('$p_{T,jet}$')
plt.xticks(np.arange(-0.5, 6.5), [400, 450, 500, 550, 600, 650, 700])
plt.tight_layout()
# plt.savefig('../../../Plots/peak_position_comparison.png', format='png', dpi=300)
plt.show()
