import numpy as np
import matplotlib.pyplot as plt

N = 150        
steps = 300     
T = 1.2         
chain = np.random.choice([1, -1], size=N)
history_1d = np.zeros((steps, N))

for t in range(steps):
    for _ in range(N):
        i = np.random.randint(0, N)
        neighbors = chain[(i-1)%N] + chain[(i+1)%N]
        dE = 2 * chain[i] * neighbors
        
        if dE <= 0 or np.random.rand() < np.exp(-dE / T):
            chain[i] *= -1
            
    history_1d[t] = chain.copy()

grid = np.random.choice([1, -1], size=(N, N))
history_2d = np.zeros((steps, N))

for t in range(steps):
    for parity in [0, 1]:
        i_idx, j_idx = np.indices((N, N))
        mask = (i_idx + j_idx) % 2 == parity
        
        top = np.roll(grid, 1, axis=0)
        bottom = np.roll(grid, -1, axis=0)
        left = np.roll(grid, 1, axis=1)
        right = np.roll(grid, -1, axis=1)
        
        dE = 2 * grid * (top + bottom + left + right)
        flip = (dE <= 0) | (np.random.rand(N, N) < np.exp(-dE / T))
        
        grid[mask & flip] *= -1
        
    history_2d[t] = grid[N//2, :].copy()

fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=150)
fig.suptitle(f"Space-Time Diagrams at Low Temperature ($T={T}$)", 
             fontsize=16, fontweight='bold')

# 1D Panel
axs[0].imshow(history_1d, cmap='binary', aspect='auto', interpolation='nearest')
axs[0].set_title("1D System (The Ising Chain)", fontweight='bold')
axs[0].set_xlabel("Position along the chain ($x$)")
axs[0].set_ylabel("Time (Monte Carlo Steps)")

# 2D Panel
axs[1].imshow(history_2d, cmap='binary', aspect='auto', interpolation='nearest')
axs[1].set_title("2D System (Onsager Transverse Slice)", fontweight='bold')
axs[1].set_xlabel("Position along the transverse row ($x$)")
axs[1].set_ylabel("Time (Monte Carlo Steps)")

fig.text(
    0.27, 0.02,
    "1D: Thermal noise wins.\n"
    "Domain walls wander randomly\n"
    "and no long-range order emerges.",
    ha="center",
    fontsize=10,
    style='italic',
    bbox={"facecolor": "white", "alpha": 0.9, "pad": 5}
)

fig.text(
    0.75, 0.02,
    "2D: Topology wins.\n"
    "Orthogonal connections anchor the spins,\n"
    "freezing the system into macroscopic order.",
    ha="center",
    fontsize=10,
    style='italic',
    bbox={"facecolor": "white", "alpha": 0.9, "pad": 5}
)

plt.tight_layout(rect=[0, 0.08, 1, 0.95])  
plt.savefig("comparison_1D_vs_2D.png", bbox_inches='tight')
plt.show()