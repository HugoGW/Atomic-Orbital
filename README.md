# Atomic-Orbital of an hydrogen atom
In this code, we plot the atomic structure of a hydrogen atom which is given by the wave function $\psi$.

I resolved analyticaly the Schrödinger equation with the central potential : $V \equiv V(r)$ which is :

$$\displaystyle i \hbar \frac{\partial \psi}{\partial t} = - \frac{\hbar^2}{2m} \Delta \psi - \frac{e^2}{4\pi^2\varepsilon_0 r}$$ 

and the solution is given by :

$\psi \equiv \psi_{n,l,m}(r,\theta,\phi) = R_{n,l}(r)Y_{l,m}(\theta,\phi)$ where $R_{n,l}(r)= \sqrt{ \big(\frac{2Z}{na_0} \big)^3 \frac{(n-l-1)!}{2n(n+l)!}} e^{-\frac{Zr}{na_0}} \big(\frac{2Z}{na_0} \big)^l L_{n-l-1}^{2l+1} \big(\frac{2Z}{na_0} \big)$ is the radial function of $\psi_{n,l,m}$ and $Y_{l,m}(\theta,\phi) = (-1)^m \sqrt{\frac{(2l+1)!(l-m)!}{4\pi(l+m)!}} P_l^m\big(cos(\theta) \big) e^{im\phi}$ is the orbital function.

where $n \in \mathbb{N}, l<n$ and $\lvert m \rvert \leq l$

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




