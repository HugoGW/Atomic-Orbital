import numpy as np
from numpy.linalg import norm
from scipy.special import eval_genlaguerre as laguerrel, sph_harm as spherical_harmonics
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import RegularGridInterpolator
from matplotlib.animation import FuncAnimation

Lx = Ly = Lz = 30

Nx, Ny, Nz = 50, 50, 50

a0 = 1

n, l, m = 4, 3, 0

def ψ(x, y, z, n, l, m):
    def R(r, n, l):
        # Radial wavefunction for hydrogen atom.
        ρ = 2 * r/a0 / n
        return np.exp(-ρ) * ρ**l * laguerrel(n-1-l, 2*l+1, 2*ρ)

    def Y(l, m, θ, φ):
        # Spherical harmonic function for hydrogen atom.
        return abs(spherical_harmonics(m, l, θ, φ))

    r = norm((x, y, z))
    θ = np.arctan2(np.sqrt(x**2 + y**2), z)
    φ = np.arctan2(y, x)

    return 4 * np.pi * r**2 * R(r, n, l)**2 * Y(l, m, θ, φ)**2


def orbital_2D_pic(n, l, m):
    '''
    Plots a 2D image of the orbital defined by nlm.
    '''
    Psi = np.zeros((Nx, Ny))
    X = np.linspace(-Lx/2, Lx/2, Nx)
    Y = np.linspace(-Ly/2, Ly/2, Ny)
    for i, x in enumerate(X):
        for j, y in enumerate(Y):
            Psi[i, j] = ψ(x, y, 0, n, l, m)
        
    Psi /= np.trapz(np.trapz(Psi,X),Y)
    
    plt.imshow(Psi, cmap='inferno', extent=[-Lx/2, Lx/2, -Ly/2, Ly/2], interpolation='nearest', norm=Normalize(vmin=np.min(Psi), vmax=np.max(Psi)))
    plt.text(-Lx/2.1, -Ly/2.1, f'n={n} l={l} m={m}', fontsize=13, color='white')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar(label='Probability density')
    plt.show()


def orbital_2D_movie(n, l, m):
    z = -Lz/2
    dz = Lz/Nz
    Psi = np.zeros((Nx, Ny))
    X = np.linspace(-Lx/2, Lx/2, Nx)
    Y = np.linspace(-Ly/2, Ly/2, Ny)

    fig, ax = plt.subplots(figsize=(10, 8))
    pcm = ax.pcolormesh(X, Y, Psi, cmap='inferno', shading='auto', vmin=0, vmax=1)
    plt.colorbar(pcm, ax=ax, label='Probability density')
    plt.text(-Lx/2.1, -Ly/2.1, f'n={n} l={l} m={m}', fontsize=13, color='white')
    plt.xlabel('x')
    plt.ylabel('y')

    while z < Lz/2:
        for i, x in enumerate(X):
            for j, y in enumerate(Y):
                Psi[i, j] = ψ(x, y, z, n, l, m)
        
        Psi /= np.trapz(np.trapz(Psi, X), Y)
        
        pcm.set_array(Psi)
        pcm.set_clim(vmin=np.min(Psi), vmax=np.max(Psi))
        plt.pause(0.0001)
        
        z += dz
    
    plt.show()



def orbital_3D_pic(n, l, m, num_points=10000):
    '''
    Plots a 3D image of the orbital defined by nlm.
    '''
    X = np.linspace(-Lx/2, Lx/2, Nx)
    Y = np.linspace(-Ly/2, Ly/2, Ny)
    Z = np.linspace(-Lz/2, Lz/2, Nz)
    
    Psi = np.zeros((Nx, Ny, Nz))
    for i, x in enumerate(X):
        for j, y in enumerate(Y):
            for k, z in enumerate(Z):
                Psi[i, j, k] = ψ(x, y, z, n, l, m)
    
    # Normalize the wavefunction
    Psi /= np.trapz(np.trapz(np.trapz(Psi, X), Y), Z)
    
    # Flatten the arrays for random sampling
    X_flat, Y_flat, Z_flat = np.meshgrid(X, Y, Z, indexing='ij')
    X_flat = X_flat.flatten()
    Y_flat = Y_flat.flatten()
    Z_flat = Z_flat.flatten()
    Psi_flat = Psi.flatten()
    
    # Filter points with x > 0
    mask = Z_flat > 0
    X_flat = X_flat[mask]
    Y_flat = Y_flat[mask]
    Z_flat = Z_flat[mask]
    Psi_flat = Psi_flat[mask]
    
    # Random sampling according to probability density Psi
    probabilities = Psi_flat / Psi_flat.sum()
    sampled_indices = np.random.choice(len(X_flat), size=num_points, p=probabilities)
    
    X_sampled = X_flat[sampled_indices]
    Y_sampled = Y_flat[sampled_indices]
    Z_sampled = Z_flat[sampled_indices]
    Psi_sampled = Psi_flat[sampled_indices]
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the points
    sc = ax.scatter(X_sampled, Y_sampled, Z_sampled, c=Psi_sampled, cmap='inferno', marker='.', s=5)
    
    ax.set_box_aspect([2, 2, 1])
    
    plt.title(f'Orbital: n={n}, l={l}, m={m}')
    plt.xlabel('x')
    plt.ylabel('y')
    ax.set_zlabel('z')
    plt.colorbar(sc, ax=ax, label='Probability density')
    plt.show()




orbital_2D_pic(n, l, m)
#orbital_2D_movie(n, l, m)
orbital_3D_pic(n, l, m)

