import numpy as np
import matplotlib.pyplot as plt

simXest = np.loadtxt("../build/simXest.txt")
simX = np.loadtxt("../solver/mhe/simX.txt")

x0_bar = np.array([0.0, np.pi, 0.0, 0.0])

NX = 4
y_label = [r"$x$", r"$\theta$", r"$v$", r"$\dot{\theta}$"]

print('difference |x0_est - x0_bar|', np.linalg.norm(x0_bar - simXest[0, :]))
print('difference |x_est - x_true|', np.linalg.norm(simXest - simX))


for i in range(NX):
    plt.subplot(NX,1,i+1)
    plt.plot(simX[:, i], label="true", marker="o", fillstyle="none")
    plt.plot(simXest[:, i], label="estimated", marker="+")
    plt.legend()
    plt.ylabel(y_label[i])
    plt.grid()

plt.show()
