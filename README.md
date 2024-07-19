# Atomic-Orbital of an hydrogen atom
In this code, we plot the atomic structure of a hydrogen atom which is given by the wave function $\psi$.

I resolved analyticaly the Schrödinger equation with the central potential : $V \equiv V(r)$ which is :

$i \hbar \frac{\partial \psi}{\partial t} = - \frac{\hbar^2}{2m} \Delta \psi - \frac{e^2}{4\pi^2\varepsilon_0 r}$ and the solution is given by :

$\psi \equiv \psi_{n,l,m} = R_{n,l}(r)Y_{l,m}(\theta,\phi)$ where $R_{n,l}(r)= \sqrt{ \bigg(\frac{2Z}{na_0} \bigg)^3 \frac{(n-l-1)!}{2n(n+l)!}} e^{-\frac{Zr}{na_0}} \bigg(\frac{2Z}{na_0} \bigg)^l L^{2l+1}_{n-l-1} \bigg(\frac{2Z}{na_0} \bigg)$ is the radial function of $\psi_{n,l,m}$ and $Y_{l,m}(\theta,\phi)$ is the orbital function.


