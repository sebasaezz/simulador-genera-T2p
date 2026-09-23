# Idea base
Se va a crear un simulador general que esta vez va a basar la señal en un vector de magnetización. 

Para un `Isochromat` se va a almacenar la magnetización transversal
$$
M_{xy}=M_{x}+jM_{y},
$$
así como la del eje $z$
$$
M_{z}.
$$
Además, se almacena la variable fija de off-resonance microscópico $\Delta \omega$.
![[esfericas.png|352]]

La fase transversal será
$$
\phi(t) = \operatorname{arg}\left( M_{xy}(t) \right)  \equiv -\Delta \omega t\quad (\operatorname{mod}2\pi ) .
$$
La señal a simular será la magnetización transversal
$$
S(t)= \frac{1}{N} \sum_{i=1}^{N}M_{xy,\,i}(t)= 
\frac{1}{N} \sum_{i=1}^{N}\left\lvert M_{xy,\,i}(t) \right\rvert \mathrm{e}^{ -\mathrm{i} \Delta \omega t }.
$$
El ángulo polar se va a definir como $\theta$.

El punto inicial de la magnetización será
$$
\mathbf{M}=\begin{bmatrix}
0 \\
0 \\
M_{0}
\end{bmatrix},
$$
luego se le deberán poder aplicar un `HardPulse` que simule vagamente ser un pulso RF de [[Nutación en MRI|flip angle]] $\alpha$ que haga
$$
\mathbf{M}'=R_{x}(\alpha)\mathbf{M},
$$
con $R_{x}$ la [[Matriz de rotación por un ángulo|matriz de rotación]] en el eje $x$.

Para la implementación de [[takeshima.pdf|Takeshima]], se deberá guardar
$$
\partial_{\omega}\mathbf{M},
$$
pero eso se pospone.


# Simulación convencional

Se logra hacer una secuencia se echo simple.

![[echo convencional.png]]