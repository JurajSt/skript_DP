
a = 45000.0/60/60




import matplotlib.pyplot as plt
plt.plot([2,2,3,4], [1,2,3,4], 'ro')
plt.ylabel('some numbers')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def f(t):
    return np.exp(-t) * np.cos(2*np.pi*t)

t1 = np.arange(0.0, 5.0, 0.1)
t2 = np.arange(0.0, 5.0, 0.02)

plt.figure(2)
plt.subplot(210)
plt.plot(t1, f(t1), t2, f(t2), 'k')

plt.subplot(212)
plt.plot(t2, np.cos(2*np.pi*t2), 'r--')
plt.savefig('graf.pdf')