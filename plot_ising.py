import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
N = 256                  
frames = 300            
T_start = 3.0           
T_end = 0.5            
Tc = 2.269              
sweeps_per_frame = 5    

T_vals = np.linspace(T_start, T_end, frames)
lattice = np.random.choice([1, -1], size=(N, N))


def mcmc_step(mat, T):
    for _ in range(N * N):
        i = np.random.randint(0, N)
        j = np.random.randint(0, N)
        s = mat[i, j]
        nb = mat[(i+1)%N, j] + mat[(i-1)%N, j] + mat[i, (j+1)%N] + mat[i, (j-1)%N]
        dE = 2 * s * nb
        if dE <= 0 or np.random.rand() < np.exp(-dE / T):
            mat[i, j] = -s

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5), gridspec_kw={'width_ratios': [4, 1]})
fig.suptitle("Ising model simulation", fontsize=16, fontweight='bold')


im = ax1.imshow(lattice, cmap='binary', vmin=-1, vmax=1, interpolation='nearest')
ax1.set_title("Spin structure")
ax1.axis('off')


ax2.set_xlim(0, 1)
ax2.set_ylim(T_end, T_start)
ax2.set_xticks([]) X
ax2.set_ylabel("Temperature ($T$)", fontsize=12)


ax2.axhline(Tc, color='red', linestyle='--', linewidth=2)
ax2.text(1.1, Tc, f'Tc={Tc}', color='red', va='center', fontweight='bold')


rect = plt.Rectangle((0.2, T_end), 0.6, T_start - T_end, color='royalblue', alpha=0.7)
ax2.add_patch(rect)

plt.tight_layout()


def update(frame):
    T_actual = T_vals[frame]
    

    for _ in range(sweeps_per_frame):
        mcmc_step(lattice, T_actual)
    

    im.set_data(lattice)
    

    rect.set_height(T_actual - T_end)
    

    if T_actual < Tc:
        rect.set_color('firebrick')
    else:
        rect.set_color('royalblue')

    if frame % 15 == 0:
        print(f"Progreso: {int(frame/frames*100)}% | T = {T_actual:.2f}")
        
    return im, rect
ani = animation.FuncAnimation(fig, update, frames=frames, interval=50, blit=True)

ani.save('ising_termometro.gif', writer='pillow', fps=20,dpi=250)
