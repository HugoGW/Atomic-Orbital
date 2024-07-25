# Atomic-Orbital of an hydrogen atom
In this code, we plot the atomic structure of a hydrogen atom which is given by the wave function $\psi$.

I solved analytically the Schrödinger equation with the central potential : $V \equiv V(r)$ which is :

$$\displaystyle i \hbar \frac{\partial \psi}{\partial t} = - \frac{\hbar^2}{2m} \Delta \psi - \frac{e^2}{4\pi^2\varepsilon_0 r}$$ 

and the solution is given by :

$\psi \equiv \psi_{n,l,m}(r,\theta,\phi) = R_{n,l}(r)Y_{l,m}(\theta,\phi)$ where $R_{n,l}(r)= \sqrt{ \big(\frac{2Z}{na_0} \big)^3 \frac{(n-l-1)!}{2n(n+l)!}} e^{-\frac{Zr}{na_0}} \big(\frac{2Z}{na_0} \big)^l L_{n-l-1}^{2l+1} \big(\frac{2Z}{na_0} \big)$ is the radial function of $\psi_{n,l,m}$ and 

$Y_{l,m}(\theta,\phi) = (-1)^m \sqrt{\frac{(2l+1)!(l-m)!}{4\pi(l+m)!}} P_l^m\big(cos(\theta) \big) e^{im\phi}$ is the orbital function.

where $n \in \mathbb{N}^*$, $l \lt n$ and $\lvert m \rvert \leq l$

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

$\textbf{I - Picture of an orbital in 2D}$

In this part, we're going to use the function ψ(x, y, z, n, l, m) written above to plot a picture of dimension $N_x \times N_y$ of an atomic orbital of an hydrogenoid.

First, we define our wave function $\textit{Psi}$ as a matrix of zeros of the size $N_x \times N_y$ and we create an array of $X$ and $Y$ :

    Psi = np.zeros((Nx, Ny))
        X = np.linspace(-Lx/2, Lx/2, Nx)
        Y = np.linspace(-Ly/2, Ly/2, Ny)

For each i and j of X and Y, we plot a dot of a certain color that matches the probability of an electron being at the coordinate (i,j), then we normalize it and plot the result.

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


$\textbf{II - Video of an orbital in 2D}$

Now, we're going to use and modify the code in the $\textbf{Part I}$ above. In fact, in the "loop" we wrote $\textit{Psi[i, j] = ψ(x, y, 0, n, l, m)}$ which means that we plot $ψ(x, y, 0, n, l, m)$ for $z=0$. Thus, from now on, we're going to plot $ψ(x, y, z, n, l, m)$ $\forall z \in [-L_z/2, L_z/2]$ with a space length $N_z$ and a space step $dz = L_z/N_z$.

    
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

So, how does the loop work ?

We define our wave function $\textit{Psi}$ like before and we define $z=-L_z/2$ that will be used in the loop : While $z \lt L_z/2$ $\longrightarrow$ we plot the wave function $\textit{Psi[i, j] = ψ(x, y, z, n, l, m)}$, then we make a pause and we plot the wave function for $z+dz$ : $\textit{Psi[i, j] = ψ(x, y, z+dz, n, l, m)}$ etc. Until $z + N_z dz=L_z/2$

$\textbf{III - Picture of an orbital in 3D}$

The function orbital_3D_pic(n, l, m, num_points=10000) plots a 3D image of the orbital defined by the quantum numbers $n$, $l$, and $m$.

We begin by creating 3 arrays $X, Y$ and $Z$ to represent the spatial dimensions, ranging from $-L/2$ to $L/2$ and the matrix $\textit{Psi}$ with dimensions $N_x\times N_y \times N_z$:

    X = np.linspace(-Lx/2, Lx/2, Nx)
        Y = np.linspace(-Ly/2, Ly/2, Ny)
        Z = np.linspace(-Lz/2, Lz/2, Nz)
        
        Psi = np.zeros((Nx, Ny, Nz))

For each dot $(i,j,k) \in (x,y,z)$ we plot a colored dot that matches with the probability of an electron being at this coordinate, then we normalize by using the trapezoidal integration method in all three spatial dimensions.

In order to see the inside of the 3D atomic orbital, a mask is applied to filter out points where $z \leq 0$. This is done using the condition Z_flat > 0.

The flattened and filtered probability densities are normalized to create a probability distribution, then random sampling is performed using $\textit{np.random.choice} with the probabilities provided by the normalized probability densities. This ensures that points with higher probability densities are more likely to be sampled.

The indices of the sampled points are used to extract the corresponding $X$, $Y$, $Z$, and $\textit{Psi}$ values.

To sum up, we create a new 3D grid of points for $X$, $Y$ and $Z$ with 

    X_flat, Y_flat, Z_flat = np.meshgrid(X, Y, Z, indexing='ij')

Then, we convert the 3D grid arrays into 1D arrays for easier manipulation and sampling :

    X_flat = X_flat.flatten()
    Y_flat = Y_flat.flatten()
    Z_flat = Z_flat.flatten()
    Psi_flat = Psi.flatten()

The mask "Z_flat > 0" creates a boolean array where only points with $Z>0$ are True and applying this mask to X_flat, Y_flat, Z_flat, and Psi_flat filters out the points where $z \leq 0$.

The probability density array $\textit{Psi}$_$\textit{flat}$ is normalized by dividing it by its sum, creating a valid probability distribution where the sum of all probabilities is 1.

Then, $\textit{np.random.choice}$ uses the probability distribution to sample indices from the filtered arrays. This ensures that points with higher probability density have a higher chance of being chosen.

Finally, the sampled points (X_sampled, Y_sampled, Z_sampled, Psi_sampled) are used to create a 3D scatter plot and the color of each point in the scatter plot corresponds to its probability density, providing a visual representation of regions with higher electron probability.

$\textbf{IV - Results}$

I lightly modified the code to gather all the plots and animations in 1 window with animation with the 3D plot to see how does the orbital look like.



https://github.com/user-attachments/assets/3ccc2d5e-1316-4bbe-986b-99481d5744ed


https://github.com/user-attachments/assets/5807affc-0b5c-4786-9c66-5f1442cb775b


https://github.com/user-attachments/assets/f7b871aa-0af3-43bf-9d50-43c46247cc04


https://github.com/user-attachments/assets/2c9bd6b2-c357-45dc-923b-822d47123291


PS : I'm aware that my code isn't optimized, even though I have some ideas for improvement. Therefore, feel free to modify it to reduce the calculation time.




