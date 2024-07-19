# Atomic-Orbital of an hydrogen atom
In this code, we plot the atomic structure of a hydrogen atom which is given by the wave function $\psi$.

I resolved analyticaly the Schrödinger equation with the central potential : $V \equiv V(r)$ which is :

$$\displaystyle i \hbar \frac{\partial \psi}{\partial t} = - \frac{\hbar^2}{2m} \Delta \psi - \frac{e^2}{4\pi^2\varepsilon_0 r}$$ 

and the solution is given by :

$\psi \equiv \psi_{n,l,m}(r,\theta,\phi) = R_{n,l}(r)Y_{l,m}(\theta,\phi)$ where $R_{n,l}(r)= \sqrt{ \big(\frac{2Z}{na_0} \big)^3 \frac{(n-l-1)!}{2n(n+l)!}} e^{-\frac{Zr}{na_0}} \big(\frac{2Z}{na_0} \big)^l L_{n-l-1}^{2l+1} \big(\frac{2Z}{na_0} \big)$ is the radial function of $\psi_{n,l,m}$ and $Y_{l,m}(\theta,\phi) = (-1)^m \sqrt{\frac{(2l+1)!(l-m)!}{4\pi(l+m)!}} P_l^m\big(cos(\theta) \big) e^{im\phi}$ is the orbital function.

where $n \in {1,2...}$, $l<n$ and $\lvert m \rvert \leq l$

This part is given in the code by : 

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

$\textbf{I - Orbital in 2D}$

In this part, we're going to use the function ψ(x, y, z, n, l, m) written above to plot a picture of dimension $N_x \times N_y$ of an atomic orbital of an hydrogenoid.

First, we define our wave function Psi as a matrix of zeros of the size $N_x \times N_y$ and we create an array of $X$ and $Y$ :

    Psi = np.zeros((Nx, Ny))
        X = np.linspace(-Lx/2, Lx/2, Nx)
        Y = np.linspace(-Ly/2, Ly/2, Ny)

For each i and j of X and Y, we plot a dot of a certain color that matches with the probability of an electron to be at the coordinate (i,j) then we normalize it and plot the result.

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


![image](https://github.com/user-attachments/assets/97a31dee-6771-4ee2-a055-86f2cb71dbf5)








