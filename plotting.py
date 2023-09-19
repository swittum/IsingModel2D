import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import h5py

print('Creating animation')

with h5py.File('./lattice.h5', 'r') as file:
    animation = np.array(file['animation'])
    magnetization = np.array(file['magnetization'])
    energy = np.array(file['energy'])

    fig, ax = plt.subplots(figsize=(5, 5))

    def animate(i):
        ax.clear()
        ax.axis('off')
        ax.imshow(animation[50*i])
        return []
    

    ani = FuncAnimation(fig, animate, frames=round(len(animation)/50), interval=10, blit=True)
    writer = FFMpegWriter(fps=30, metadata=dict(artist='Me'), bitrate=1800)

    # Save the animation as a video
    ani.save('./ani.mp4', writer=writer)
    plt.close(fig)
    
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))
    for ax in (ax1, ax2):
        ax.grid()
    ax1.set_title('Magnetization')
    ax1.plot(magnetization)
    ax2.set_title('Energy')
    ax2.plot(energy)
    plt.tight_layout()
    plt.savefig('./magnetization_energy.png')
    # plt.show(block=True)
        
    

    

    
