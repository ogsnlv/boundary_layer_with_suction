# boundary_layer_with_suction
The code describes the self-similar solution of 2D boundary-layer equations for a flow around an infinite plate using tridiagonal matrix method.

The govergning equations are  


$$
u\frac{\partial u}{\partial x} + v\frac{\partial u}{\partial y}=\nu \frac{\partial^2 u}{\partial^2 y}
$$,
$$
\frac{\partial u}{\partial x} + \frac{\partial v}{\partial y}=0;
$$


with boundary conditions:

$$u = 0, \ v=V_0(x) \ \text{at} \ y=0;
$$ u \rightarrow U_\infty \ \text{at} \  y \rightarrow \infty$$

where $V_0(x)$ represents the suction/injection at the plate surface. 

The calculations are performed in C#, whereas the plots are created with Python and **matplotlib**.


