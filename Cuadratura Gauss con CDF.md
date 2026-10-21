La integral es
$$
\int_{-\infty}^{\infty} m(t,\omega)p(\omega) \, \mathrm{d}\omega 
$$
tal que
$$
p(\omega) = \frac{a}{\pi[\omega^{2}+(1/a)^{2}]},\quad a=1/T_{2}'.
$$
Con CDF
$$
F(\omega)= \frac{1}{2}+\frac{1}{\pi}\arctan\left( \frac{\omega}{a} \right)\implies \omega = F^{-1}(u)=a\tan \left[ \pi \left( u-\frac{1}{2} \right)  \right] .
$$
Se hace el cambio:
$$
S(t)=\int_{0}^{1} m(t, F^{-1}(u)) \, \mathrm{d}u 
$$
La integral queda como
$$
S(t)=\int_{0}^{1} m(t, \textstyle a\tan \left[ \pi \left( u-\frac{1}{2} \right)  \right] ) \, \mathrm{d}u
$$
Usando el cambio
$$
x = 2u-1\implies \mathrm{d}u = \frac{1}{2}\mathrm{d}x
$$
permitiendo cambiar los límites a
$$
S(t)=\frac{1}{2}\int_{-1}^{1} m\left( t, a \tan\left( \frac{\pi x}{2} \right) \right) \, \mathrm{d}x 
$$

Esto permite usar cuadratura Gauss-Legendre
$$
S(t)=\frac{1}{2}\sum_{i=0}^{N} v_{i} \cdot  m\left( t, a \tan\left( \frac{\pi x_{i}}{2} \right) \right).
$$
