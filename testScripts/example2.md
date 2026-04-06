<!-- Name: Mesh generation in CFD: All you need to know in one article -->
<!-- Source: https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/mesh-generation-in-cfd-all-you-need-to-know-in-one-article/ -->

## In this series

[custom_category_posts_list category_slug="10-key-concepts-everyone-must-understand-in-cfd"]

## In this article

[TOC]

## Introduction

I have a confession to make. For the longest time, I have avoided linear systems of equations like the pest. Even when I had to solve them, like the Pressure Poisson equation, I used whatever method was quickest to implement, without giving much thoughts towards convergence properties or performance. I was very happy to wait for a month for the simulations to finish, and I was living a happy life without matrices in my life.

I think this has to do with my (professional) upbringing; as a student, you are very impressionable, and one thing I often heard is that linear system of equations solvers, i.e. those that solve $\mathbf{Ax}=\mathbf{b}$, are, well, linear. The governing equations of fluid mechanics are non-linear, and so, the reason was that, as we have to linearise the system of equations, we are losing some non-linearity.

Well, this is true, but, we are able to recover the non-linearity through an iterative procedure (something which we will touch upon at the end of this article), and while it will cost additional iterations, the convergence rate that you will obtain from solving $\mathbf{Ax}=\mathbf{b}$ usually offsets this additional cost.

In this article, I want to give you an introduction to methods to solve $\mathbf{Ax}=\mathbf{b}$ that I was always lacking. The combination of fearmongering by some of my professors and the sense of comfort explicit methods gave me meant that I never really felt compelled to look into implicit systems in earnest. But once I did, I regretted having waited for so long.

In this article, I want to equip you with all the necessary knowledge so that you can integrate implicit solvers, i.e. those that solve $\mathbf{Ax}=\mathbf{b}$, into your CFD codes.

Before we jump into all of the various methods to solve $\mathbf{Ax}=\mathbf{b}$, let us first look at how this linear system of equations comes about, and then we will look at methods to solve it!

### A unified view on how Ax=b is constructed and solved

In the following, I want to give you a unified view on how $\mathbf{Ax}=\mathbf{b}$ is constructed. I haven't seen this written anywhere (perhaps it is and I just haven't seen it), but after I saw linear systems from this unified perspective, I felt that I had a much better and deeper, intuitive understanding of how to construct explicit and implicit time integration scheme.

To start us off, we need to talk about [explicit](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/explicit-vs-implicit-time-integration-and-the-cfl-condition/#aioseo-explicit-time-integration) and [implicit](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/explicit-vs-implicit-time-integration-and-the-cfl-condition/#aioseo-implicit-time-integration) time integration. An explicit time integration scheme is one where, after discretisation, only one unknown is left in the discretised equation. We can solve for this unknown *explicitly*. Conversely, an implicit time integration scheme is one where, after discretisation, more than one unknowns are obtained, and we are only able to solve this through a linear system of equations, i.e. by solving $\mathbf{Ax}=\mathbf{b}$.

Let's take a simple, 1D equation to see what I mean. Let us discretise the following unsteady heat-diffusion equation:

$$
\frac{\partial T}{\partial t}=\Gamma\frac{\partial^2 T}{\partial x^2}
\tag{eq:unsteady-heat-diffusion}
$$

Let's discretise that with a simple first-order Euler method in time, and a second-order central scheme in space. We will denote the unknown (next) time level by $n+1$, while all variables at the current time level $n$ are known. Therefore, we obtain the following equations:

- **Explicit** time integration:

$$
\frac{T_i^{n+1}-T_i^{n}}{\Delta t} = \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}
\tag{eq:explicit-heat-diffusion}
$$

- **Implicit** time integration:

$$
\frac{T_i^{n+1}-T_i^{n}}{\Delta t} = \Gamma\frac{T^{n+1}_{i+1} - 2T^{n+1}_{i} + T^{n+1}_{i-1}}{(\Delta x)^2}
$$

We see that the only difference is at which time level we evaluate the temperatures on the right-hand side, i.e. either at $n$ or at $n+1$. The next step is to collect all unknowns on the left-hand side of the equation, all knowns on the right-hand side of the equation, and then to write the equation in coefficient form. This is done in the following way for both schemes:

- **Explicit** time integration:

$$
\frac{T_i^{n+1}}{\Delta t} = \frac{T_i^{n}}{\Delta t} + \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}\\[1em]
\left[\frac{1}{\Delta t}\right]T_i^{n+1} = \frac{T_i^{n}}{\Delta t} + \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}
\tag{eq:explicit-system-discretised}
$$

- **Implicit** time integration:

$$
\frac{T_i^{n+1}}{\Delta t} - \Gamma\frac{T^{n+1}_{i+1} - 2T^{n+1}_{i} + T^{n+1}_{i-1}}{(\Delta x)^2} = \frac{T_i^{n}}{\Delta t}\\[1em]
\left[\frac{1}{\Delta t}\right]T_i^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i+1}^{n+1} + \left[\frac{2\Gamma}{(\Delta x)^2}\right]T_{i}^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i-1}^{n+1} = \frac{T_i^{n}}{\Delta t}\\[1em]
\left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i+1}^{n+1} + \left[\frac{1}{\Delta t} + \frac{2\Gamma}{(\Delta x)^2}\right]T_{i}^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i-1}^{n+1} = \frac{T_i^{n}}{\Delta t}
\tag{eq:implicit-system-discretised}
$$

A common notation here is to introduce coefficients $a_P$, $a_E$, and $a_W$ for the coefficients in-front of $T^{n+1}_{i}$, $T^{n+1}_{i+1}$, and $T^{n+1}_{i-1}$, respectively, with $E$ and $W$ standing for east and west. Furthermore, all of the known quantities on the right-hand side can be put together into a single value, which we will denote as $b$. Doing so gives us the following system of equations:

$$
a_E \cdot T^{n+1}_{i+1} + a_P \cdot T^{n+1}_{i} + a_W \cdot T^{n+1}_{i-1} = b
\tag{eq:coefficient-form-for-lionear-system}
$$

We can then write the coefficients as:

- **Explicit** time integration:

$$
a_E = 0,\qquad a_P=\frac{1}{\Delta t},\qquad a_W = 0,\qquad b = \frac{T_i^{n}}{\Delta t} + \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}
$$

- **Implicit** time integration:

$$
a_E = \frac{-\Gamma}{(\Delta x)^2}, \qquad a_P=\frac{1}{\Delta t} + \frac{2\Gamma}{(\Delta x)^2}, \qquad a_W = \frac{-\Gamma}{(\Delta x)^2}, \qquad b = \frac{T_i^{n}}{\Delta t}
$$

Since this Eq.(\ref{eq:coefficient-form-for-lionear-system}) is given in terms of $i$, we will obtain this equation for each node $i$ in our mesh. So, if we have $nx=5$ nodes in our 1D mesh, then we will obtain Eq.(\ref{eq:coefficient-form-for-lionear-system}) for 5 different nodes.

We have to take special care at boundaries, that is, how we impose values on the left and right side of the domain. We can either impose Dirichlet, Neumann, or Robin boundary conditions, and how to do that, including for linear systems like the ones we are currently looking at, are subject of the next article on [boundary conditions](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-implement-boundary-conditions-in-cfd/).

Let's say we want to impose Dirichlet boundary conditions on the left and right side of the domain; that is, we want to impose a fixed temperature on the left and right side of the domain. Let's say we set the temperature to $T_L=0$ on the left boundary of the domain, and $T_R=1$ on the right boundary of the domain. If we want to do that, the trick here is to set $a_P=1$ and $a_W=a_E=0$, as well as $b_L=T_L$ and $b_R=T_R$. To make this clear, on the left side of the domain, at the boundary, we impose:

$$
a_E = 0,\qquad a_P=1,\qquad a_W = 0,\qquad b = T_L
$$

On the right side of the domain, we have:

$$
a_E = 0,\qquad a_P=1,\qquad a_W = 0,\qquad b = T_R
$$

If we insert these values back into Eq.(\ref{eq:coefficient-form-for-lionear-system}), we obtain (now using $T_{L,R}$ to write the equation for both left and right sides of the domain):

$$
0 \cdot T^{n+1}_{i+1} + 1 \cdot T^{n+1}_{i} + 0 \cdot T^{n+1}_{i-1} = T_{L,R}\\[1em]
T^{n+1}_{i} = T_{L,R}
$$

I have gone over this a bit quicker, but I have a much more in-depth discussion for how to apply these types of boundary conditions in the next article in the section [How to implement boundary conditions for implicit time integration](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-implement-boundary-conditions-in-cfd/#aioseo-how-to-implement-boundary-conditions-for-implicit-time-integration). Give that a read if you want to have a more thorough explanation.

With these boundary conditions defined, we can now write all of the equations, for each node $i$, in matrix form, that is, we can write what we have just discretised as $\mathbf{Ax}=\mathbf{b}$. Let's stay with the example of using 5 nodes, this means we have 2 nodes on the boundary, and 3 nodes on the internal domain. A sketch of this domain is shown in the following:

<!-- wp:image {"width":"500px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\1D_domain.png" alt="Sketch of a 1D domain with boundary conditions on the left and right, as well as internal nodes" class="wp-image-5550" style="width:500px"/></figure>
<!-- /wp:image -->

For the 1D domain, noting that $T=\mathbf{x}$, we can now write our linear system $\mathbf{Ax}=\mathbf{b}$ as:

$$
\begin{bmatrix}
a_P & a_E & 0 & 0 & 0 \\[1em]
a_W & a_P & a_E & 0 & 0 \\[1em]
0 & a_W & a_P & a_E & 0 \\[1em]
0 & 0 & a_W & a_P & a_E \\[1em]
0 & 0 & 0 & a_W & a_P \\[1em]
\end{bmatrix}
\begin{bmatrix}
x_0 \\[1em]
x_1 \\[1em]
x_2 \\[1em]
x_3 \\[1em]
x_4 \\[1em]
\end{bmatrix}=
\begin{bmatrix}
b_0 \\[1em]
b_1 \\[1em]
b_2 \\[1em]
b_3 \\[1em]
b_4 \\[1em]
\end{bmatrix}
\tag{eq:linear-system-general}
$$

We can now insert values we have for our coefficients $a_P$, $a_E$, and $a_W$, as well as for the right-hand side vector $b$. For the first and last row, we have to apply our boundary conditions, that is, $a_P=1$, $a_E=a_W=0$, and $b_{L,R}=T_{L,R}$. This is true for both explicit and implicit time integrations. Let's do this now. For both explicit and implicit time integrations, we obtain the following linear systems of equations:

- **Explicit** time integration:

$$
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
0 & 1/\Delta t & 0 & 0 & 0 \\[1em]
0 & 0 & 1/\Delta t & 0 & 0 \\[1em]
0 & 0 & 0 & 1/\Delta t & 0 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
\begin{bmatrix}
x_0 \\[1em]
x_1 \\[1em]
x_2 \\[1em]
x_3 \\[1em]
x_4 \\[1em]
\end{bmatrix}=
\begin{bmatrix}
T_L \\[1em]
\frac{T_1^{n}}{\Delta t} + \Gamma\frac{T^{n}_{2} - 2T^{n}_{1} + T^{n}_{0}}{(\Delta x)^2} \\[1em]
\frac{T_2^{n}}{\Delta t} + \Gamma\frac{T^{n}_{3} - 2T^{n}_{2} + T^{n}_{1}}{(\Delta x)^2} \\[1em]
\frac{T_3^{n}}{\Delta t} + \Gamma\frac{T^{n}_{4} - 2T^{n}_{3} + T^{n}_{2}}{(\Delta x)^2} \\[1em]
T_R \\[1em]
\end{bmatrix}
\tag{eq:linear-system-explicit-example}
$$

- **Implicit** time integration:

$$
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
-\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 & 0 & 0 \\[1em]
0 & -\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 & 0 \\[1em]
0 & 0 & -\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
\begin{bmatrix}
x_0 \\[1em]
x_1 \\[1em]
x_2 \\[1em]
x_3 \\[1em]
x_4 \\[1em]
\end{bmatrix}=
\begin{bmatrix}
T_L \\[1em]
\frac{T_1^{n}}{\Delta t} \\[1em]
\frac{T_2^{n}}{\Delta t} \\[1em]
\frac{T_3^{n}}{\Delta t} \\[1em]
T_R \\[1em]
\end{bmatrix}
\tag{eq:linear-system-implicit-example}
$$

To see that this really does represent our original equations, you can take any row in the coefficient matrix, multiply that by the solution vector $\mathbf{x}$, which represents the temperature in our example, and set that equal to the right-hand side vector $\mathbf{b}$. For example, taking the third row in both cases results in:

- **Explicit** time integration:

$$
\left[\frac{1}{\Delta t}\right]x_2^{n+1} = \frac{x_2^{n}}{\Delta t} + \Gamma\frac{T^{n}_{1} - 2T^{n}_{2} + T^{n}_{3}}{(\Delta x)^2}\\[1em]
\left[\frac{1}{\Delta t}\right]x_i^{n+1} = \frac{T_i^{n}}{\Delta t} + \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}
$$

- **Implicit** time integration:

$$
\left[\frac{-\Gamma}{(\Delta x)^2}\right]x_{1}^{n+1} + \left[\frac{1}{\Delta t} + \frac{2\Gamma}{(\Delta x)^2}\right]x_{2}^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]x_{3}^{n+1} = \frac{T_2^{n}}{\Delta t}\\[1em]
\left[\frac{-\Gamma}{(\Delta x)^2}\right]x_{i+1}^{n+1} + \left[\frac{1}{\Delta t} + \frac{2\Gamma}{(\Delta x)^2}\right]x_{i}^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]x_{i-1}^{n+1} = \frac{T_i^{n}}{\Delta t}
$$

I have given both versions here, i.e. with explicit indices $i=1,2,3$, and with $i-1, i, i+1$. If we can now make the mental substitution of $x_{1,2,3}=T_{1,2,3}$, then we see that the above equations really are the same as Eq.(\ref{eq:explicit-system-discretised}) and Eq.(\tag{eq:implicit-system-discretised}).

OK, so now we need to solve our linear system of equations. From linear algebra, we know that we can solve $\mathbf{Ax}=\mathbf{b}$ by inverting the coefficient matrix and then left-multiplying it with the equation, i.e.:

$$
\mathbf{Ax}=\mathbf{b}\\[1em]
\mathbf{A}^{-1}\mathbf{Ax}=\mathbf{A}^{-1}\mathbf{b}\\[1em]
\mathbf{Ix}=\mathbf{A}^{-1}\mathbf{b}\\[1em]
\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}
$$

Here, $\mathbf{I}$ is the identity matrix, i.e. containing ones on the main diagonal and zeros everywhere else. Let's look at the explicit time integration case first. Our coefficient matrix $\mathbf{A}$ is:

$$
\mathbf{A}=
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
0 & 1/\Delta t & 0 & 0 & 0 \\[1em]
0 & 0 & 1/\Delta t & 0 & 0 \\[1em]
0 & 0 & 0 & 1/\Delta t & 0 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
$$

Since we only have elements on the diagonal, we can easily invert this matrix with the following rule:

$$
a'_{ii}=1/a_{ii}
$$

Here, $a'_{ii}$ are the entries in the inverted matrix, i.e. $\mathbf{A}^{-1}$ at row $i$ and column $i$. Since both row and column indices are always the same, we will always pick the value on the diagonal. The values of $a_{ii}$ are coming from the diagonal of the original coefficient matrix $\mathbf{A}$. With this rule, we can now invert the coefficient matrix as:

$$
\mathbf{A}^{-1}=
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
0 & \Delta t & 0 & 0 & 0 \\[1em]
0 & 0 & \Delta t & 0 & 0 \\[1em]
0 & 0 & 0 & \Delta t & 0 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
$$

We can also write our discretised system for $\mathbf{Ax}=\mathbf{b}$ now as $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$, this becomes:

$$
\begin{bmatrix}
x_0 \\[1em]
x_1 \\[1em]
x_2 \\[1em]
x_3 \\[1em]
x_4 \\[1em]
\end{bmatrix}=
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
0 & \Delta t & 0 & 0 & 0 \\[1em]
0 & 0 & \Delta t & 0 & 0 \\[1em]
0 & 0 & 0 & \Delta t & 0 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
\begin{bmatrix}
T_L \\[1em]
\frac{T_1^{n}}{\Delta t} + \Gamma\frac{T^{n}_{2} - 2T^{n}_{1} + T^{n}_{0}}{(\Delta x)^2} \\[1em]
\frac{T_2^{n}}{\Delta t} + \Gamma\frac{T^{n}_{3} - 2T^{n}_{2} + T^{n}_{1}}{(\Delta x)^2} \\[1em]
\frac{T_3^{n}}{\Delta t} + \Gamma\frac{T^{n}_{4} - 2T^{n}_{3} + T^{n}_{2}}{(\Delta x)^2} \\[1em]
T_R \\[1em]
\end{bmatrix}
\tag{eq:linear-system-explicit-example}
$$

Again, taking the third row as an example, we can write out the discretised equation as:

$$
x_2^{n+1}=\Delta t\left[\frac{T_2^{n}}{\Delta t} + \Gamma\frac{T^{n}_{3} - 2T^{n}_{2} + T^{n}_{1}}{(\Delta x)^2}\right]\\[1em]
x_2^{n+1}=\Delta t\frac{T_2^{n}}{\Delta t} + \Delta t\Gamma\frac{T^{n}_{3} - 2T^{n}_{2} + T^{n}_{1}}{(\Delta x)^2}\\[1em]
x_2^{n+1}=T_2^{n} + \frac{\Gamma\Delta t}{(\Delta x)^2}\left[T^{n}_{3} - 2T^{n}_{2} + T^{n}_{1}\right]\\[1em]
x_i^{n+1}=T_i^{n} + \frac{\Gamma\Delta t}{(\Delta x)^2}\left[T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}\right]
\tag{eq:explicit-heat-diffusion-for-node-two}
$$

Again, if we can mentally replace $x_2=T_2$ in this case, then we realise that this is really just Eq.(\ref{eq:explicit-heat-diffusion}), which we saw at the very beginning. We can easily rearrange Eq.(\ref{eq:explicit-heat-diffusion}) to look like Eq.(\ref{eq:explicit-heat-diffusion-for-node-two}). So, all of this discussion of bringing Eq.(\ref{eq:explicit-heat-diffusion}) into the form of $\mathbf{Ax}=\mathbf{b}$, inverting the coefficient matrix, and then solving for $\mathbf{x}$ seems a bit over the top, doesn't it?

Hold on to that thought, we will be coming back to that in just a second! But, before we do, let's also look at the implicit time integration case quickly. Here, the coefficient matrix was given as:

$$
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
-\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 & 0 & 0 \\[1em]
0 & -\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 & 0 \\[1em]
0 & 0 & -\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
$$

In this case, we have elements both on the diagonal and on the off-diagonal. Because of these off-diagonal entries in $\mathbf{A}$, we can no longer find the inverse, i.e. $\mathbf{A}^{-1}$, in the same way as we did with the explicit time integration case. Thus, we are forced to compute it properly, and, as it turns out, once we have large systems, this becomes impossible.

Just to bring home this point: The number of nodes or cells in your mesh determines the size of your coefficient matrix. If you have 1 million cells in your mesh, your matrix will have 1 million rows and columns. Inverting a 1 million by 1 million matrix, even if most of the entries are zero, is essentially impossible, and we have to potentially do that at each timestep and/or iteration.

To add insult to injury, 1 million cells, in the context of CFD, is a pretty small simulation. You may have some simple test cases that have only a few hundred thousand cells. For example, a classical 2D airfoil simulation will probably have about 100,000 cells. But any 3D engineering applications that are of industrial interest can easily reach the hundreds of millions of cells. In one case, I am aware of a company routinely running 1 billion cells per simulation. 

If you are planning to write a CFD solver for these types of applications and you propose to invert the matrix at each iteration, I think you may be shot. Don't know, but that feels about right. Don't do it.

There are now essentially 2 things we can do: either we try to numerically approximate the inverse of $\mathbf{A}^{-1}$, or we just try to find a corrective algorithm that will try, starting from some initial guess $\mathbf{x}^(0)$, to find a solution that satisfies our original system of equations $\mathbf{Ax}=\mathbf{b}$. Typically, we do the latter.

The rest of this article is essentially dedicated to exploring these methods, which find the solution of $\mathbf{x}$ in an iterative procedure. These iterative procedures mainly differ in their approach to solving the linear system, but some are better than others in terms of convergence and stability.

Thankfully, we have a lot of experience with all of these different methods, and while the literature is filled with all sorts of different algorithms, there really are just a few that we need to know about. These are the ones you will find in actual CFD solvers. All of the other methods may have some relevance outside the field of CFD, but we do not need to take a great interest in them.

Before we start looking at some matrix properties, I wanted to come back to my unified view, which I was telling you about before. I showed you that both explicit and implicit time integrations can be represented as $\mathbf{Ax}=\mathbf{b}$. Typically, when we write an explicit solver, we would probably not write the explicit system in this form, simply because it takes more lines of code and complexity.

But there is one very important reason why I would encourage you to even write your system like that. If you do, you can use any matrix vector library that will do all of the computation for you. That, on its own, may not be necessarily an advantage, but if that library offers parallelisation support, you can use that without having to think about parallelising the code yourself.

And, if you do decide to support implicit time stepping at some point in the future, well, the infrastructure will already be there in your code; you just have to rederive the system of equations for that particular discretisation. Well, even that can be automated by packages like Sympy, Maple, or even Matlab's symbolic toolbox.

Hopefully, this makes sense. So, if you see $\mathbf{Ax}=\mathbf{b}$, just know that this really represents the governing equation you are trying to solve. Change the coefficient matrix $\mathbf{A}$ or the right-hand side vector $\mathbf{b}$, and you solve a different equation. With that out of the way, let's start looking at some matrix properties.

## Matrix properties

At this point in the discussion, we have to talk about matrices. If you just think of them as a 2D storage container, oh boy, you are in for an awakening. People have filled entire books about matrices alone. Believe me, people have made a lot of effort trying to study and demistify matrices. So, in this section, we will look at some of the their properties.

### Dense and sparse matrices

Let's start off easy and look at the two structures a matrix can be in. A dense matrix is one in which all, or almost all entries in the matrix are non-zero. In contrast, a sparse matrix contains mostly zeros, where typically only the diagonal and some off-diagonals contain non-zero entries. The off-diagonals are shifted from the diagonal of the matrix either to the right or down. The following sketch shows this, where a cross indicates that this entry is a non-zero value.

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\dense_vs_sparse_matrix.png" alt="Comparison between a dense and a sparse matrix, shopwing that a dense matrix has all (or almost all) entries filled in the matrix, while a sparse matrix will have only diagonal and off-diagonal contributions" class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

If we go back to our example in the previous section, we had the following matrix $\mathbf{A}$ for the implicit time integration:

$$
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
-\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 & 0 & 0 \\[1em]
0 & -\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 & 0 \\[1em]
0 & 0 & -\Gamma/(\Delta x)^2 & 1/\Delta t + 2\Gamma/(\Delta x)^2 & -\Gamma/(\Delta x)^2 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
$$

This matrix contains elements only on the diagonal (either $1$ or $1/\Delta t + 2\Gamma/(\Delta x)^2$), while the two off-diagonals contain the entries $-\Gamma/(\Delta x)^2$. So, this is a sparse matrix, and all of our matrices that arise in CFD (as part of discretising our governing equations) will be sparse. Sparse matrices themselves are a topic you could do an entire PhD on (if you were *that* committed), and we will touch upon some key properties in a second.

### Lower and Upper triangular matrices

In the analysis that will follow, we will make heavy use of lower (L) and upper (U) triangular matrices. If you already know, or have heard, of the LU decomposition, that is what the letter stands for.

Since we can add matrices together, that means that we can decompose them into as many matrices we want. A common decomposition is to split a matrix into its lower, upper, and sometimes also diagonal part, as shown in the following sketch for a dense matrix:

<!-- wp:image {"width":"800px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\LUD_decomposition.png" alt="This sketch shows how a matrix can be decomposed into a lower triangular matrix L, a diagonal matrix D (containing only elements on the diagonal), and an upper triangular matrix U, containing only elements on the upper left triangle of the matrix." class="wp-image-5550" style="width:800px"/></figure>
<!-- /wp:image -->

We can thus say that if we are dealing with a matrix $\mathbf{A}$, that we can write this matrix as well as:

$$
\mathbf{A} = \mathbf{L} + \mathbf{D} + \mathbf{U}
$$

Sometimes, we write this decomposition simply as:

$$
\mathbf{A} = \mathbf{L} + \mathbf{U}
$$

If we do that, then the diagonal matrix has been either absorbed into the upper or lower triangular matrix.

### Symmetric and non-symmetric matrices

Now that we can decompose our matrix into its lower and upper triangular matrix, we can easily check if our matrix is symmetric or not. A symmetric matrix is one where all the elemens are mirrored along the diagonal of a matrix. That means that if I have a matrix entry $a_{2,4}$, the coefficient in a matrix at row 2 and column 4, it will be the same as $a_{4,2}$, i.e. the coefficient at row 4 and column 2. We can generalise this as $a_{i,j}=a_{j,i}$.

Since the transpose of a matrix is essentially just flipping a matrix at its diagonal, we can say that the transpose of a symmetric matrix is just the matrix itself. In mathematical terms, we can write:

$$
\mathbf{A}^T=\mathbf{A}
$$

And, since we have introduced the decomposition into the lower and upper triangular part already, we can also look at what happens if we transpose either the lower or upper triangular matrix. In this case, I refer to the lower and upper triangular matrices which do not contain the diagonal, i.e. the $\mathbf{D}$ is treated separately.

If we take the lower triangular matrix and transpose it, we flip it at the diagonal. So it will now occupy the space where the upper triangular matrix is. If the matrix is symemtric, then both of these will have the same values, i.e. we transpose of the lower triangular matrix is the same as the upper triangular matrix. The same is true for the transpose of the upper triangular matrix, which is the same as the lower triangular matrix. We can write this as:

$$
\mathbf{L}^T=\mathbf{U}\\
\mathbf{U}^T=\mathbf{L}
$$

### How matrices are stored

For a dense matrix, we do have to store the entire matrix, i.e. each coefficient $a_{i,j}$. There is no way around that. For small problem sizes, that isn't a problem, but consider that a matrix will have $m\cdot n$ entries. The matrices that we will look at will all be square matrices, so $m=n$. This means that our matrix will contain $m^2$ entries. So, if I have a mesh with 10 cells, my matrix will contain $m^2=10^2=100$ entries. If I double the number of cells to 20, then my matrix will contain $m^2=20^2=400$.

Thus, the number of elements we have to store grow exponentially (to the power of 2), and if we are dealing with 1 million elements in our mesh, well, we have to store $m^2=1,000,000^2 = 1,000,000,000,000$. We don't do that, especially if our matrix is sparse. If we had a sparse matrix, which contains one diagonal and, for the sake of argument, let's say we also store 6 off-diagonal entries. Then, we store $1,000,000$ entries in the diagonal term, and a few less in each off-diagonal.

If we look back to our dense vs. sparse matrix sketch above, we saw that the diagonal in the sparse matrix contains 8 entries, while there are two off-diagonals that contain 7 entries and two off-diagonals that contain 4 entries. So, the further the off-diagonals are away from the diagonal, the fewer entries they will have, but in our example, if the off-diagonals have a a few elements missing compared to the diagonal (which has $1,000,000$ entries), then they will still have close to $1,000,000$ entries. 

Let's just assume that each off-diagonal also stores $1,000,000$. Then, we stores one diagonal + 6 off-diagonals, so $7,000,000$ entries in total. If we stored that in a dense matrix, where each zero would be stored as well, we said that we would be needing to store $1,000,000,000,000$ entries. If we take the ratio, than we can estimate how many elements in that matrix will actually be non-zero. We have:

$$
\frac{7,000,000}{1,000,000,000,000}=0.000007=0.0007\%=7\cdot 10^{-4}\%
$$

So, we only need 0.0007% of the matrix, the remaining 99.9993% are just zeros. For this reason, we really want to store sparse matrices in a different from, not as full, dense matrix.

There are a few different forms for how we can store sparse matrices, and conceptually the easiest to understand is the coordinate matrix, as shown in the following image: 

<!-- wp:image {"width":"800px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\coordinate_matrix.gif" alt="An animation showing how each element in a sparse matrix is stored within the coordinate storage format." class="wp-image-5550" style="width:800px"/><figcaption class="wp-element-caption">Figure reproduced from <a href="https://matteding.github.io/2019/04/25/sparse-matrices/" target="_blank" rel="noopener" title="">Matt Eding - Sparse Matrices</a>.</figcaption></figure>
<!-- /wp:image -->

We need to store the row, the column, and the value at this location. We do that by storing all of this information in three separate vectors. If we are now evaluation a matrix vector product of the form $\mathbf{Ab}$, where $\mathbf{A}$ is the matrix and $\mathbf{b}$ being a vector, then we will have to multiply each row of $\mathbf{A}$ with the full vector $\mathbf{b}$. If we say that we count rows with index $i$ and columns with index $j$, then we can say that each row will perform the following computation:

$$
\sum_{j=0}^{nColumns}a_{i,j}b_j
\tag{eq:matrix-vector-multiplication-for-one-row}
$$

So, if we wanted to evaluate this with a sparse matrix which we store in coordinate form, we would have to check, for each row, if any values are actually stored. So, refering to the figure above, if we are in row 5, we see that there are no entries here. There is also no 5 in our row vector, as shown to the right of the matrix. Thus, we can say that, since there is no data here, the result of the summation in Eq.(\ref{eq:matrix-vector-multiplication-for-one-row}) is zero.

However, if we go to row 3, for example, then we do have a value at column location $j=4$. Since we start with a zero index, $j=4$ corresponds to the fifth column, and so we have to multiply this by the fifth entry in the vector $\mathbf{b}$, i.e. $b_4$.

So, when we evaluate Eq.(\ref{eq:matrix-vector-multiplication-for-one-row}), we have to check now, for each value $a_{i,j}$, if it exists. If it does, great, we can perform thje computation. If it doesn't, then we set $a_{i,j}=0$. In reality, we would implement the summation as shown in Eq.(\ref{eq:matrix-vector-multiplication-for-one-row}), though, as this would be very wasteful.

Let's go back to our previous example of the 1 million by 1 million coefficient matrix. Here we said that we had $1,000,000,000,000$ entries. So, we would not be checking $1,000,000,000,000$ times if an entry existed, that would be just as bad as storing the matrix itself. To make this clear, let's look at the example of evaluating $\mathbf{c}=\mathbf{Ab}$.

To make this a concrete example, let's use some numbers. I will be using the same matrix as shown above in the figure, but I have come up with a very creative vector, if I may say myself, that we'll use to multiply our matrix with (i.e. the vector $\mathbf{b}$). The resulting vector $\mathbf{c}$ can then be calculated as:

<!-- wp:image {"width":"400px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\coordinate_matrix_example.png" alt="An example calculation of " class="wp-image-5550" style="width:400px"/></figure>
<!-- /wp:image -->

Look at Eq.(\ref{eq:matrix-vector-multiplication-for-one-row}) again, for each row in the vector $\mathbf{c}$, we have to sum over all of the products of $\mathbf{A}$ and $\mathbf{b}$, where the coefficients in $\mathbf{A}$ are non-zero. For the first row, we see that the thrid column has a non-zero value of 9, so we multiply that by the third entry in $\mathbf{b}$, which is 3, and so the corresponding product is $9\cdot 3 = 27$. We add that to $\mathbf{c}$ at the first row.

Since we do not have any additional coefficient in $\mathbf{A}$ in the first row, we do not add anything else to $\mathbf{c}$ in the first row. Well, this is just matrix vector multiplication, hopefully this is something you have heard before. And if not, you must be a child under the age of 18.

I recently had to update my privacy policy, and I used one of those privacy policy generator, that asks you a few questions and then it will customise the privacy policy for you. I looked through it, and according to that privacy policy, children (which were counting as those that are under the age of 18), are not allowed to *consume* the content on this website.

Well, what harm is there in a 17 year old, a day away from their 18th birthday, to learn about matrix vector multiplication I ask? I didn't quite see the point. I removed that paragraph and now anyone from the age of 0-99+ is allowed to *consume* my content to their heart's content. You are welcome.

In any case, now that we have an example, let's implement that into code. I have written a small python script that will do this for us. Have a look through it, I'll explain it below if it is not clear.

<!-- wp:kevinbatdorf/code-block-pro {"code":"# coordinate matrix form\nrows = [1, 3, 0, 2, 4]\ncolumns = [1, 4, 2, 3, 3]\ndata = [2, 5, 9, 1, 6]\n\nnumber_of_non_zero_elements = len(data)\n\n# vector b which we use in matrix vector multiplication\nb = [1, 2, 3, 4, 5, 6]\n\n# result vector c which is initialised to 0 everywhere\nc = [0, 0, 0, 0, 0, 0]\n\nfor i in range(0, number_of_non_zero_elements):\n    row_id = rows[i]\n    column_id = columns[i]\n    c[row_id] += data[i] * b[column_id]\n\nprint(c) # this prints: [27, 4, 4, 25, 24, 0]","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# coordinate matrix form\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003erows = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e4\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003ecolumns = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e4\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003edata = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e9\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e6\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003enumber_of_non_zero_elements = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003elen\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(data)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# vector b which we use in matrix vector multiplication\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eb = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e4\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e6\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# result vector c which is initialised to 0 everywhere\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003ec = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003ein\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003erange\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, number_of_non_zero_elements):\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    row_id = rows\u0026#91;i\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    column_id = columns\u0026#91;i\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    c\u0026#91;row_id\u0026#93; += data\u0026#91;i\u0026#93; * b\u0026#91;column_id\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eprint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(c) \u003c/span\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# this prints: \u0026#91;27, 4, 4, 25, 24, 0\u0026#93;\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"python","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":19,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(2 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly># coordinate matrix form
rows = &#91;1, 3, 0, 2, 4&#93;
columns = &#91;1, 4, 2, 3, 3&#93;
data = &#91;2, 5, 9, 1, 6&#93;

number_of_non_zero_elements = len(data)

# vector b which we use in matrix vector multiplication
b = &#91;1, 2, 3, 4, 5, 6&#93;

# result vector c which is initialised to 0 everywhere
c = &#91;0, 0, 0, 0, 0, 0&#93;

for i in range(0, number_of_non_zero_elements):
    row_id = rows&#91;i&#93;
    column_id = columns&#91;i&#93;
    c&#91;row_id&#93; += data&#91;i&#93; * b&#91;column_id&#93;

print(c) # this prints: &#91;27, 4, 4, 25, 24, 0&#93;</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #6A9955"># coordinate matrix form</span></span>
<span class="line"><span style="color: #D4D4D4">rows = &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">4</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">columns = &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">4</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">data = &#91;</span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">9</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">6</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">number_of_non_zero_elements = </span><span style="color: #DCDCAA">len</span><span style="color: #D4D4D4">(data)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># vector b which we use in matrix vector multiplication</span></span>
<span class="line"><span style="color: #D4D4D4">b = &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">4</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">6</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># result vector c which is initialised to 0 everywhere</span></span>
<span class="line"><span style="color: #D4D4D4">c = &#91;</span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"></span>
<span class="line"><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> i </span><span style="color: #C586C0">in</span><span style="color: #D4D4D4"> </span><span style="color: #DCDCAA">range</span><span style="color: #D4D4D4">(</span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, number_of_non_zero_elements):</span></span>
<span class="line"><span style="color: #D4D4D4">    row_id = rows&#91;i&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">    column_id = columns&#91;i&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">    c&#91;row_id&#93; += data&#91;i&#93; * b&#91;column_id&#93;</span></span>
<span class="line"></span>
<span class="line"><span style="color: #DCDCAA">print</span><span style="color: #D4D4D4">(c) </span><span style="color: #6A9955"># this prints: &#91;27, 4, 4, 25, 24, 0&#93;</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

First, we are going to set all values in $\mathbf{c}$ to zero. This is important, as we will be adding results to this vector now, and we may be adding more than once to the same location in $\mathbf{c}$. Then, we will loop over the number of elements in the coefficient matrix $\mathbf{A}$. If we look at the example calculation above, we had 5 entries in total in $\mathbf{A}$, so our loop length would be 5 (```number_of_non_zero_elements```), and not $6\cdot 7=42$, which is the size of the matrix shown (with 6 rows and 7 columns).

Then at each loop iteration, we get the current row, the current column, and the value at that location. We now also get the value in our $\mathbf{b}$ vector, specifically at the same location as our column location that we got from our matrix. We multiply these values together, and then add to $\mathbf{c}$, specifically at the same row location we received from the matrix. We need to add this value to whatever was previously stored in this vector, since we may have more than one value per row stored in $\mathbf{A}$.

As you can see, the results that are being printed are the ones we would expect. Now, this is a quick and dirty example, by just hardcoding the rows, columns, and data entries, but we could also make this a bit cleaner by creating a ```CoordinateMatrix``` class, which stores these arrays. If we then overload the multiplication operator, i.e. the ```*``` operator, then we could write our matrix vector product as ```c = A * b```. Would't that be nice? Well, that isn't really a big issue in Python, and we could achieve this by writing:

<!-- wp:kevinbatdorf/code-block-pro {"code":"class CoordinateMatrix:\n    def __init__(self):\n        self.rows = []\n        self.columns = []\n        self.data = []\n    \n    def __mul__(self, b):\n        c = [0] * len(b)\n        for i in range(0, len(self.rows)):\n            row_id = self.rows[i]\n            column_id = self.columns[i]\n            c[row_id] += self.data[i] * b[column_id]\n        return c\n\n\nA = CoordinateMatrix()\nA.rows = [1, 3, 0, 2, 4]\nA.columns = [1, 4, 2, 3, 3]\nA.data = [2, 5, 9, 1, 6]\n\nb = [1, 2, 3, 4, 5, 6]\nc = A * b\n\nprint(c) # prints [27, 4, 4, 25, 24, 0]","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eclass\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e \u003c/span\u003e\u003cspan style=\u0022color: #4EC9B0\u0022\u003eCoordinateMatrix\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e:\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003edef\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003e__init__\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e):\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e.rows = []\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e.columns = []\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e.data = []\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003edef\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003e__mul__\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003eb\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e):\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        c = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93; * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003elen\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(b)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003ein\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003erange\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003elen\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e.rows)):\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e            row_id = \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e.rows\u0026#91;i\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e            column_id = \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e.columns\u0026#91;i\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e            c\u0026#91;row_id\u0026#93; += \u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eself\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e.data\u0026#91;i\u0026#93; * b\u0026#91;column_id\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003ereturn\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e c\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eA = CoordinateMatrix()\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eA.rows = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e4\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eA.columns = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e4\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eA.data = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e9\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e6\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eb = \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e3\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e4\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e6\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003ec = A * b\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eprint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(c) \u003c/span\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# prints \u0026#91;27, 4, 4, 25, 24, 0\u0026#93;\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"python","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":24,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(2 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly>class CoordinateMatrix:
    def __init__(self):
        self.rows = []
        self.columns = []
        self.data = []
    
    def __mul__(self, b):
        c = &#91;0&#93; * len(b)
        for i in range(0, len(self.rows)):
            row_id = self.rows&#91;i&#93;
            column_id = self.columns&#91;i&#93;
            c&#91;row_id&#93; += self.data&#91;i&#93; * b&#91;column_id&#93;
        return c


A = CoordinateMatrix()
A.rows = &#91;1, 3, 0, 2, 4&#93;
A.columns = &#91;1, 4, 2, 3, 3&#93;
A.data = &#91;2, 5, 9, 1, 6&#93;

b = &#91;1, 2, 3, 4, 5, 6&#93;
c = A * b

print(c) # prints &#91;27, 4, 4, 25, 24, 0&#93;</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #569CD6">class</span><span style="color: #D4D4D4"> </span><span style="color: #4EC9B0">CoordinateMatrix</span><span style="color: #D4D4D4">:</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #569CD6">def</span><span style="color: #D4D4D4"> </span><span style="color: #DCDCAA">__init__</span><span style="color: #D4D4D4">(</span><span style="color: #9CDCFE">self</span><span style="color: #D4D4D4">):</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #569CD6">self</span><span style="color: #D4D4D4">.rows = []</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #569CD6">self</span><span style="color: #D4D4D4">.columns = []</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #569CD6">self</span><span style="color: #D4D4D4">.data = []</span></span>
<span class="line"><span style="color: #D4D4D4">    </span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #569CD6">def</span><span style="color: #D4D4D4"> </span><span style="color: #DCDCAA">__mul__</span><span style="color: #D4D4D4">(</span><span style="color: #9CDCFE">self</span><span style="color: #D4D4D4">, </span><span style="color: #9CDCFE">b</span><span style="color: #D4D4D4">):</span></span>
<span class="line"><span style="color: #D4D4D4">        c = &#91;</span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">&#93; * </span><span style="color: #DCDCAA">len</span><span style="color: #D4D4D4">(b)</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> i </span><span style="color: #C586C0">in</span><span style="color: #D4D4D4"> </span><span style="color: #DCDCAA">range</span><span style="color: #D4D4D4">(</span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #DCDCAA">len</span><span style="color: #D4D4D4">(</span><span style="color: #569CD6">self</span><span style="color: #D4D4D4">.rows)):</span></span>
<span class="line"><span style="color: #D4D4D4">            row_id = </span><span style="color: #569CD6">self</span><span style="color: #D4D4D4">.rows&#91;i&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">            column_id = </span><span style="color: #569CD6">self</span><span style="color: #D4D4D4">.columns&#91;i&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">            c&#91;row_id&#93; += </span><span style="color: #569CD6">self</span><span style="color: #D4D4D4">.data&#91;i&#93; * b&#91;column_id&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #C586C0">return</span><span style="color: #D4D4D4"> c</span></span>
<span class="line"></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">A = CoordinateMatrix()</span></span>
<span class="line"><span style="color: #D4D4D4">A.rows = &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">4</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">A.columns = &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">4</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">A.data = &#91;</span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">9</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">6</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">b = &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">3</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">4</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">6</span><span style="color: #D4D4D4">&#93;</span></span>
<span class="line"><span style="color: #D4D4D4">c = A * b</span></span>
<span class="line"></span>
<span class="line"><span style="color: #DCDCAA">print</span><span style="color: #D4D4D4">(c) </span><span style="color: #6A9955"># prints &#91;27, 4, 4, 25, 24, 0&#93;</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

If this makes sense, then congratulation, you have mastered the most difficult bit already about sparse matrices and how we store them. However, the coordinate form is not as efficient as some other sparse matrix storage formats, and so we don't actually use this format in practice, at least not for the actual matrix vector multiplication.

The format you will most likely use when you are dealing with sparse matrices is the compressed sparse row, or CSR format. Now, I did already talk about the CSR format in depth in my article on [sparse matrices in CFD application](https://cfd.university/learn/how-to-compile-write-and-use-cfd-libraries-in-c/how-to-write-a-cfd-library-the-sparse-matrix-class/#aioseo-sparse-matrices-in-cfd-applications), which you may want to consult. In that article I also show how we can implement this format into a C++ code. But I want to give a brief overview here as well, so that you understand the basic concept.

The CSR format is conceptually very similar to the coordinate matrix form (sometimes also abbreviated as the COO form). The following figure shows how we store a matrix in CSR form:

<!-- wp:image {"width":"800px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\compressed_row_storage.gif" alt="An animation showing how each element in a sparse matrix is stored within the compressed sparse row storage format." class="wp-image-5550" style="width:800px"/><figcaption class="wp-element-caption">Figure reproduced from <a href="https://matteding.github.io/2019/04/25/sparse-matrices/" target="_blank" rel="noopener" title="">Matt Eding - Sparse Matrices</a>.</figcaption></figure>
<!-- /wp:image -->

We have again three arrays that we are storing, but we store different information here. Arguably, the information that we do store is organised in a way that isn't as clear as the coordinate form, and there are good reasons for it. So, let's look at this storage form and then discuss why this is better.

First up, we have the ```index pointer``` array, and it has the same length as the number of rows + 1. We will get to this array in a second. We also have an ```indices``` array, which will give us the column indices of where values are store, and then we have a ```data``` array again, which stores the actual non-zero values.

OK, so we have a ```indices``` array, which will give us the column index, but where is the array giving us the corresponding row entry? Well, it is encoded in the ```index pointer``` array. The index pointer array, we said, had a length of the number of rows + 1. We can use this fact and say that if we want to get some information about, say, row 0, or 4, all I have to do is to go into the corresponding location within the ```index pointer``` array.

OK, so let's do that. At location 0, within the ```index pointer``` array, we have a value of 0. We also need to get the next entry in this array (and for that reason, this array has one additional entry compared to the number of rows), which is 2. The location in the ```index pointer``` array that we access, in our case, that was the first location (i.e. index 0), that corresponds to the first row (again, with index 0). So what ever information we are reading, will be relevant for row 0.

We read the value at location 0, plus the next entry, and we got the values 0 and 2 from the ```index pointer``` array. The difference between the two is $2-0=2$. This tells us that there are 2 non-zero elements for row 0. We can see that this is the case from the figure above. So, we already know that we are in row 0, now we need to get the column locations (or indices). This is where the ```indices``` array comes in.

The ```index pointer``` does not only tell us how many entries there are in the corresponding row, but it will also give us the indices we have to use to get the column locations from the ```indices``` array. Which indices did we read? Well, we know that there are 2 non-zero elements in row 0, and, the ```index pointer``` at location 0 is 0 as well, so it tells us that we have to go into the ```indices``` array at location 0 (this is what we read from the ```index pointer``` array), and then we have to read 2 entries from this array.

These two entries in the ```indices``` array are 0 and 2, so we know that in row 0, there are two non-zero elements at colum 0 and 2. We can again confirm that this is correct from the figure above. We do the same trick to retrieve the non-zero values from the ```data``` array, i.e. we go to location 0 and read the next two entries here, which are 8 and 2. So, we know that at row 0 and column 0, we have a value of 8, and at row 0 and colum 2, we have a value of 2.

Let's do that again for location 4 (i.e. the fifth entry) in the ```index pointer``` array. So, for row 4, we get ```index pointer[4]=3```. The next entry in the same array gives us ```index pointer[5]=6```. This means that there are ```6-3=3``` entries in row 4. Thus, we have to go into the indices array, at l;ocation 3, and then read the next three values, which will give us the corresponding column indices/locations.

If we do that, we get ```indices[3]=2```, ```indices[4]=3```, and ```indices[5]=4```. The corresponding non-zero values are ```data[3]=7```, ```data[4]=1```, and ```data[5]=2```. Looking at the figure above, we can again see that this is correct.

Well, if we look at this format, we may be asking, why is this so much better than the coordinate matrix format? Well, for starters, we no longer store the row index explicitly. So we have already reduced the storage requirements somewhat. I'd argue that this is not really that much of a problem, for the applications that we have in mind (CFD), the coordinate format would require, perhaps, 50% more storage. So, yes, it would be larger, but not that much. If performance is our main driving force, though, CSR is an easy optimisation.

In any case, there is a much bigger advantage for CSR. In the coordinate matrix form, we were not required to store the non-zero values in any particular order. This is very lucrative when we are building the matrix and we don't know yet how our sparse matrix will look like exactly. So we can always easily append a row and column index with an associated non-zero matrix coefficient, and it doesn't matter in which order we append that.

We saw from our coordiante matrix example before that there wasn't any particular order in which the rows were stored. For the CSR matrix, though, the storage order is important, and we have to store all values row-by-row. This means that when we get data from RAM, and load that into our CPU, all values that are stored in one row will likely be within a cache line. As a result, when we perform the matrix vector multiplication, we will have fewer cache misses and so overall faster computation of our matrix vector product.

If you want to see how to implement sparse matrices using the CSR format, you can check out the implementation, as well as some additional thoughts on the CSR format in the [sparse matrices in CFD application](https://cfd.university/learn/how-to-compile-write-and-use-cfd-libraries-in-c/how-to-write-a-cfd-library-the-sparse-matrix-class/#aioseo-sparse-matrices-in-cfd-applications) article.

Some libraries that compute sparse matrix vector products may actually allow you to first assemble your matrix as a coordinate matrix (which means you can easily add non-zero values at arbitrary locations). Once you are done, the library may then lock the matrix and transform it into a CSR format, to make it more memory efficient. So, the coordinate matrix is still useful, but for performance reasons, we prefer CSR. 

I hope you can see now why we really want to store sparse matrices in a different format. If you understand both the coordinate and compressed sparse row format, you know the most important ones. There are other formats available, but they essentially do the same thing, they may jsut store different arrays and thus retrieve data differently.

### The Krylov subspace

If you are designing websites, it is customary to select a primary colour and a secondary colour for contrast. These are the only colour you will then use throughout the entire website design. Combine that with a background colour, and you have everything that you need. Can you guess which primary colour I have picked? While this website is now much more toned-down, it used to be a lot more aggressively pink (mainly because I used wordpress and I just couldn't figure out (or couldn't be asked?!) how to make the website, well, less pink).

In my cause, I have settled for a cold, dark blue, but, you may read this in the future where I have changed my mind again and adopted a tropical colour scheme, who knows. Next time you visit a website, look out for the primary and secondary colour (there is no real secondary colour on this website, at least not in the version that was live when I wrote this).

So, if we pick any colour from the rainbow of colours, we can for a subset of colours. In my case, this website uses the subset of colours which are white, dark grey, and dark blue. But we can pick any other colours as well, there is no limitation of what we can pick, apart from the colours having to look good together.

Contrast a subset now with a subspace. A subspace is similar to a subset, but, if I combine two elements from a subspace, I do get a different element that lives in the same subspace. So, if we go back to our colour example, a subspace of colours would be any shades of a specific colour. So instead of settling, for example, only for the colour of blue, I could select shades of blue, which would form both a subset and a subspace.

If I combine two shades of blue, I get another shade of blue. I can repeat this process with any shade of blue and I will always get some form of blue back. Now, if I throw in the colour green, I will loose my subspace. Blue and green will give me a different colour when combined, it will no longer be purely blue. However, blue and green still form a subset of all the colours.

So, a subset can be any set of a larger collection, and there really aren't any rules. Anything is allowed, much like the very liberal choice of the clothings my neighbour wears. I mean, is it a towl that is just too small, or is it a skirt? I can't tell (and, to be honest, I don't want to, but do I have to live with permanently closed window blinds?). Is a towl a subset of clothing? You came here for CFD, but these are the *real questions* we must contemplate in life.

A subspace, then, is like a ubset, but with some rules, or structure. In this case, combination of elements within a subspace need to produce a new element in the subspace. Let me give you a different example, which is a bit more mathematical. Let's look a tthe following vector:

$$
\begin{bmatrix}
x\\
y\\
0
\end{bmatrix}
$$

If I asked you to produce two vectors, at random, you may give me:

$$
\mathbf{v}_1=
\begin{bmatrix}
1\\
2\\
0
\end{bmatrix},\qquad
\mathbf{v}_2=
\begin{bmatrix}
3\\
4\\
0
\end{bmatrix}
$$

Excellent choice, really putting your creativity here on display. Let's see what happens when we we add these together, i.e. we compute $\mathbf{v}_1 + \mathbf{v}_2$. We get:

$$
\begin{bmatrix}
1\\
2\\
0
\end{bmatrix}+
\begin{bmatrix}
3\\
4\\
0
\end{bmatrix}=
\begin{bmatrix}
4\\
6\\
0
\end{bmatrix}
$$

Did you notice something about these vectors? Both $\mathbf{v}_1$ and $\mathbf{v}_2$ are, first of all, subsets of all the vectors in three dimensional space. Mathematically speaking, we can say that $\mathbf{v}_1$ and $\mathbf{v}_2$ are part of $\mathbb{R}^3$ (the three dimensional space), or $\mathbf{v}_1,\mathbf{v}_2\in \mathbb{R}^3$. If we want to say that $\mathbf{v}_1$ and $\mathbf{v}_2$ are a subset of $\mathbb{R}^3$, then we would write:

$$
\mathbf{v}_1, \mathbf{v}_2 \subset \mathbb{R}^3
$$

But, $\mathbf{v}_1$ and $\mathbf{v}_2$ do not just form a subset, they also form a subspace! Both $\mathbf{v}_1$ and $\mathbf{v}_2$ have $z=0$, that is, the third component of the vector is zero. Since both vectors only contain $x$ and $y$ elements, we can say that both $\mathbf{v}_1$ and $\mathbf{v}_2$ live inside the xy-plane. And, what happened when we added $\mathbf{v}_1$ and $\mathbf{v}_2$ together? The resulting vector was also part of the xy-plane.

Thus, $\mathbf{v}_1$ and $\mathbf{v}_2$ are not just a subset, but they are also a subspace! Let's look at them for a moment:

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\xy-plane_subspace.png" alt="A plot showing the v1 and v2 vector, and how these form a line in the xy-plane. Two additional vectors are show, one going along the x-axis, and one going along the y-axis. These are linear independent and map out the entire xy-plane." class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

We can see that both $\mathbf{v}_1$ and $\mathbf{v}_2$ are indeed in the xy-plane, and we can see that they span the area shown in orange. Does this mean $\mathbf{v}_1$ and $\mathbf{v}_2$ only form a subset of the subspace (i.e. we can only reach the parts shown in orange)? Not quite. As long as two vectors are linearly independent, we can reach any point in the xy plane!

We can see that both vectors are linearly independent in the figure. If they are not linear independent, then they would both point in the same direction. For example, a vector with coordinates $(1,1,0)$ and another vector with $(2,2,0)$ would be not linear independent. If we divide the second by 2, we get the first vector, so they are linear dependent.

Let's show that $\mathbf{v}_1$ and $\mathbf{v}_2$ are linear independent. We can create any vector in the xy plane as a linear combination of $\mathbf{v}_1$ and $\mathbf{v}_2$ as:

$$
\mathbf{v}=\alpha\mathbf{v}_1 + \beta\mathbf{v}_2=
\mathbf{v}=\alpha
\begin{bmatrix}
1\\
2\\
0
\end{bmatrix}+\beta
\begin{bmatrix}
3\\
4\\
0
\end{bmatrix}
$$

For example, with $\alpha=-3$ and $\beta=1.5$, we have:

$$
\mathbf{v}_3 = -3
\begin{bmatrix}
1\\
2\\
0
\end{bmatrix}+1.5
\begin{bmatrix}
3\\
4\\
0
\end{bmatrix}=
\begin{bmatrix}
-3\\
-6\\
0
\end{bmatrix}+
\begin{bmatrix}
4.5\\
6\\
0
\end{bmatrix}=
\begin{bmatrix}
1.5\\
0\\
0
\end{bmatrix}
$$

We can see this vector in the figure above and we can confirm that this vector is not just part of the orange area, but indeed outside it. Thus, using $\mathbf{v}_1$ and $\mathbf{v}_2$ alone, we can form the entire subspace, i.e. we can reach any point in the xy plane. We just have to find the values for $\alpha$ and $\beta$, that's all. 

Mathematically speaking, we can write this subspace as:

$$
\text{span}\{\mathbf{v}_1,\mathbf{v}_2\}=\{(x,y,0)\in\mathbb{R}^3\,|\,x,y\in\mathbb{R}\}
$$

Here, the span keyword states that $\mathbf{v}_1$ and $\mathbf{v}_2$ form the entire subspace. As we saw, we can reach any point in the subspace as a linear combination of these two vectors. The right-hand side states that any real numbers ($\mathbb{R}$) $x$ and $y$ form a vector $(x,y,0)$, which is part of the three dimensional space $\mathbb{R}^3$.

Let's step back and see why this subspace definition is important. Let's say we want to draw a map of a city. The city exist in three dimensional space, and some cities will have rather large elevation changes. So, we could draw an [isometric projection](https://en.wikipedia.org/wiki/Isometric_projection) of the city, show all of the elevation changes, as well as the city itself. But that's not what we normally do, is it?

Instead, we can say that we project the entire 3D city onto a 2D subspace and then we draw the city here. This simplifies our task to draw the city on a 2D piece of paper. So, certain tasks become easier in a subspace.

With this in mind, let us now turn to the Krylov subspace and why this is so important. Imagine the following scenario: You and your friend have reached fame beyond your wildest imagination, and you are making a trailer for your new and upcoming show. You both decide that jumping out of an airplane with a parachute would be a good idea, and then doing the promotion while falling through the sky. The twist: Your friend decided to prank you.

One by one, people jump out the airplane as planned. But then, the pilot comes up to you and says: "Your friend told me you always wanted to be a pilot, have fun!" and then he jumps out as well. You are left on your own and you have no idea how to fly a plane. Sounds a bit far-fetched? This actually happend in Germany. If you want to see what happened, see here: [Part1](https://www.youtube.com/watch?v=tovjgDNHLs4) and [Part2](https://www.youtube.com/watch?v=u6u9g0h8kaY) (you may want to use auto generated subtitles). As they say, *German sense of humor* ...

But let's put ourselves into the same position. We have no idea how to fly the airplane, but we need to learn how to fly it, and fast! What would you do? Well, we humans are incredible good at observing things and drawing conclusions. So, our natural instinct is to manipulate an unknown system and then see how this system responds. In this case, the unknown system is the aircraft, and so a sensible thing to do may be to push the controls forward and to pull them back to see how the aircraft responds.

We notice that pushing the controls forward, the aircraft starts to lower the nose and we start to loose altitute. As we pull back on the controls, the aircraft raises the nose again and we increase in altitute. So we have recorded two responses of the system based on our inputs. Let's now formalise this. Let's say that we have a state vector that stores our input. Let's call this vector $\mathbf{v}$. So, we push or pull the controls, and this will change the inputs to the vector $\mathbf{v}$.

This vector $\mathbf{v}$ stores the pitch (pulling or pushing the controls), the roll (turning the aircraft left or right), and the yaw (also turning the aircraft left and right but with the rudder pedals). From flight dynamics, we know that equations do exist which describe the motion of an airplane based on changes to its state (roll, pitch, and yaw). But, let's say we don't know that.

However, as we saw at the beginning of this article, we can write any equation as a matrix vector product of the coefficients from the equation and the independent variables, e.g. roll, pitch, and yaw. See Eq.(\ref{eq:linear-system-general}), for example, where we wrote our partial differential equation as a matrix vector product. So, we can conclude that if we have a vector containing roll, pitch, and yaw, we can multiply that by some matrix $\mathbf{A}$, which will then give us the response of the aircraft.

So, we pull back on the controls (we change $\mathbf{v}$), and, as we already discussed, this results in the aircraft pitching the nose up. Thus, $\mathbf{Av}$, that is, the response of the aircraft to changes in our inputs, is that the nose is raised and we start to climb. Now, here comes the next part, and this is critical. What happens to the aircraft when we are in this state of having a raised nose?

Well, if we release the controls, the aircraft will lower the nose by itslelf and go back into a state of equillibrium, this is just how aircrafts are designed (except for military jets). How can we express that in mathematical terms? We said that raising the nose and gaining altitude was a response from our aircraft as a change in $\mathbf{v}$, so we recorded that as $\mathbf{Av}$. The response of the aircraft to this state is, well, $\mathbf{A}(\mathbf{Av})$. So, we check the response of the aircraft to a state in which we brought it before.

We can also write this as: $\mathbf{A}(\mathbf{Av})=\mathbf{A}^2(\mathbf{v})$. What happens next? So, we have raised the nose, that was our input $\mathbf{Av}$. Then, we released the controls and the aircraft started to pitch down again. This was the aircraft's response to our input and this was $\mathbf{A}(\mathbf{Av})=\mathbf{A}^2(\mathbf{v})$. And now, the aircraft will initially pitch down a bit too much, and so we will actually start the decend. This is the aircraft's response to its own previous response to our change in control input.

We can write this as $\mathbf{A}(\mathbf{A}(\mathbf{Av}))=\mathbf{A}^3(\mathbf{v})$. It will now pitch up again, so it is the response to the previous state, and this could be captured as $\mathbf{A}^4(\mathbf{v})$, and so on. This motion is known as the phugoid, and it is shown in the following figure:

<!-- wp:image {"width":"800px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\phugoid_motion.jpeg" alt="This figure shows how pushing back on the control initially and then releasing the controls will result in a constant pitch up and pitch down response by the aircraft, which eventually will be dampend." class="wp-image-5550" style="width:800px"/><figcaption class="wp-element-caption">Figure reproduced from <a href="https://phugoidaero.com/" target="_blank" rel="noopener" title="">phugoidaero.com</a>.</figcaption></figure>
<!-- /wp:image -->

One of the perks of working at a university with its own airport (and its own set of research aircrafts) is that we do get to fly these type of manouvres, and it never gets boring! Depending on your weight and balance, you can either fly a stable or unstable phugoid (where the oscillation is dampend or amplified). In both cases (but in particular in the unstable mode), you get close to 0g, something which is best experienced, and difficult to describe.

In any case, let's take a step back, again, and see what we have done. We have played with the controls and we have gotten a sense for how the aircraft will react to our inputs ($\mathbf{Av}$), as well as to its own inputs, i.e. those states produces from our initial input ($\mathbf{A}^k\mathbf{v}$). We can do the same thing by turning the controls, pushing the rudders, changing the throttles, pressing buttons, etc. and then record what happens.

All of this will help us understand how the aircraft behaves. That is, we don't know anything about the internal working of the aircraft (we don't know the matrix $\mathbf{A}$), but we can test how the aircraft responses to our inputs (we can observe what the results are of $\mathbf{Av}$ and subsequent reactions of the aircraft $\mathbf{A}^k\mathbf{v}$).

Thus, our input $\mathbf{v}$ and all of multiplications with $\mathbf{A}$ form a subspace that represents how our aircraft responds. The matrix $\mathbf{A}$ is just a normal matrix, but not all matrices are part of the subspace in which the responses $\mathbf{Av}$ live. For example, if we allowed any matrix to be part of the subspace, that means that we may pull back on the controls, and instead of raising the nose, we now operate the throttles with the controls.

This is not an expected behaviour, but, if we did not have a subspace, i.e. if we said any matrix can be used here, then we would get random and unexpected behaviour. Earlier I said that a subspace has a certain structure, which differentiates it from a subset. The structure here is that changes to our pitch will always result in either nose up or nose down, but it will not suddenly operate the throttles, or the radio, or the windscreen wipers, etc.

Just for completeness, changes in pitch can, under certain circumstances, result in roll as well (pitch and roll are coupled in reality). If we fly near the stall speed, one wing may stall before the other, leading to a loss in lift on one side and an asymmetric lift distribution, which induces a roll. This, unfortunately, still leads to many  [accidences in general aviation](https://www.youtube.com/shorts/-_469TFPjGA) (small aircrafts flown for recreational purposes).

If we wanted to express that $\mathbf{v}$ and its product with $\mathbf{A}$ form a subspace, then we can write this as:

$$
\text{span}\{\mathbf{v},\mathbf{Av},\mathbf{A}^2\mathbf{v},\mathbf{A}^3\mathbf{v},...,\mathbf{A}^{k-1}\mathbf{v}\}
$$

Again, in plain english, this definition just states that each change to our inputs and the response of the aircraft to these inputs, has a predictable outcome. We can also show that with the same xy plane example we saw earlier. Let's define a matrix $\mathbf{A}$ as:

$$
\mathbf{A}=
\begin{bmatrix}
a & b & c \\
d & e & f \\
g & h & i
\end{bmatrix}
$$

Now we say that $g=h=i=0$. Then, we can write this matrix as:

$$
\mathbf{A}=
\begin{bmatrix}
a & b & c \\
d & e & f \\
0 & 0 & 0
\end{bmatrix}
$$

Let's create a specific matrix with concrete numbers. We can really just take any values here for $a$, $b$, $c$, $d$, $e$, and $f$, the choice won't matter. I am defining the following matrix:

$$
\mathbf{A}=
\begin{bmatrix}
3 & -2 & 0 \\
0 & 4 & 1 \\
0 & 0 & 0
\end{bmatrix}
$$

Let's now also define an arbitrary vector $\mathbf{v}$. Any vector will do here and I am defining:

$$
\mathbf{v}=
\begin{bmatrix}
3 \\
-2 \\
1.5
\end{bmatrix}
$$

Let's compute $\mathbf{A}\mathbf{v}$ first. This gives us:

$$
\begin{bmatrix}
3 & -2 & 0 \\
0 & 4 & 1 \\
0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
3 \\
-2 \\
1.5
\end{bmatrix}=
\begin{bmatrix}
13 \\
-6.5 \\
0
\end{bmatrix}
$$

This newly computed vector has indeed a value of zero for the third (z) component. Thus, application of our original vector $\mathbf{v}$ to our matrix $\mathbf{A}$ results in a vector that is within the xy-plane. We could say that one of the properties of $\mathbf{A}$ is that it will project vectors onto the xy-plane.

We can compute the product of $\mathbf{A}$ with the newly created vector, i.e. $\mathbf{A}(\mathbf{A}\mathbf{v})=\mathbf{A}^2\mathbf{v}$. This results in:

$$
\begin{bmatrix}
3 & -2 & 0 \\
0 & 4 & 1 \\
0 & 0 & 0
\end{bmatrix}\left(
\begin{bmatrix}
3 & -2 & 0 \\
0 & 4 & 1 \\
0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
3 \\
-2 \\
1.5
\end{bmatrix}\right)=
\begin{bmatrix}
3 & -2 & 0 \\
0 & 4 & 1 \\
0 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
13 \\
-6.5 \\
0
\end{bmatrix}=
\begin{bmatrix}
52 \\
-22 \\
0
\end{bmatrix}
$$

The resulting vector is still in the xy-plane, and so repeated applications of $\mathbf{A}$ to $\mathbf{v}$ will produce additional vectors that are within the xy-plane. Thus, we can say that repeated products of $\mathbf{A}$ with $\mathbf{v}$ will span the entire xy-plane, that is:

$$
\text{span}\{\mathbf{Av},\mathbf{A}^2\mathbf{v},\mathbf{A}^3\mathbf{v},...,\mathbf{A}^{k-1}\mathbf{v}\}
$$

Notice that $\mathbf{v}$ itself is not part of the subspace, so I have removed it from the span. When I say span, think of the entire xy-plane. In the figure we looked earlier, we saw that the vectors $\mathbf{v}_1$ and $\mathbf{v}_2$ *spaned* (created) the an area in the xy-plane which I had shown in orange. But, we saw that we could reach any point in the xy-plane by using a linear combination of $\mathbf{v}_1$ and $\mathbf{v}_2$.

This is what I mean by span, we map out an area (which in this case is a subset of the entire xy-plane), but, through linear combinations of these vectors, we can reach any other point on that plane. $\mathbf{v}_1$ and $\mathbf{v}_2$ alone are enough to go to any point on this xy-plane.

As it turns out, there is one subspace that is of great importance to us, which is the Krylov subspace. This is formally defined as:

$$
\mathcal{K}_r(\mathbf{A},\mathbf{v})=\text{span}\{\mathbf{v},\mathbf{Av},\mathbf{A}^2\mathbf{v},\mathbf{A}^3\mathbf{v},...,\mathbf{A}^{r-1}\mathbf{v}\}
$$

Why is this so important? Well, let's look at the inventor of this subspace. Alexei Nikolaevich Krylov was a naval engineer, and he tried to compute some characteristic behaviours of ships. He lived in a time in which ship dynamics wasn't that wel understood, and the combination of ships with the fluid dynamics of waves and wind wasn't all that well udnerstood.

So, what did Krylov do? Well, he *cheated*. He came up with a simplified form of the dynamics which still captured the essence of the problem, but reducing the problem significantly, so that he could find suitable approximations. A true engineer!

Let's think about another analogy. Some years ago, I travelled through Stansted Airport in the UK, with a fairly well-known, low-cost airline. No need to provide free advertisement here, so let'c call this airline Brian Air.

I arrived in good time, and I was hungry, so I went to the only place that was back then selling hot food: Burger King. I don't remember which burger I ordered, but my choice didn't have any influence anyways, as I only received the bun with the meat patty. No sauce, no salad, no flavour.

If you asked me "Is this a burger", I would have probably said yes, although that wouldn't have been a happy yes. So, we can say that the this burger, bereft of flavour, would be part of the subspace of burgers, as it had the essential building blocks (bun + patty), but not much more.

Bringing this back to Krylov's *cheating*: Burger King sold me a burger, and they technically gave me a burger (they offered me *something* resembling a burger from the burger subspace), so they cheated, too. They could have made an effort and prepare the burger according to their instructions, but they were taking some shortcuts.

Krylov's *cheating* was slightly different, in that he looked at the problem of ships and tried to capture the essence of the problem (much like the bun and patty captures the essence of a burger). Specifically, he was interested in stability and vibration of ships, and the problem he was trying to solve was very similar to the phugoid motion of the aircraft going up and down that we looked at.

So, Krylov was interested in studying the general system that describes rigid body motion of the form:

$$
\mathbf{M}\ddot{\mathbf{x}}+\mathbf{K}\mathbf{x}=f
\tag{eq:rigid-body-motion-general}
$$

Here, $\mathbf{M}$ is the mass matrix and $\mathbf{K}$ is the stiffness matrix. If we have external forces on the system, these can be expressed by $f$, otherwise we have $f=0$. For example, Krylov wanted to know something about ship vibrations, so he needed to study some form of eigenvalue problem. To do that, he had to use an assumption of what $\mathbf{x}$ would look like. A common approach back then, and still today, is to use the so-called *normal mode ansatz*.

This approach states that a vibration can be expressed through an exponential function of the form:

$$
\mathbf{x}=\mathbf{v}e^{i\omega t}
$$

Here, $\mathbf{v}$ is some constant vector, which tells us something about the modes, while $\omega$ is the angular frequency. We are using complex numbers here, so we have used $i$ here as well, and if you hate complex numbers, don't worry, we will get rid of them in a second. $t$ is the time.

Our goal is to insert our assumption for $\mathbf{x}$ (our normal mode ansatz) into Eq.(\ref{eq:rigid-body-motion-general}). To do that, we need to find the second derivative of $\mathbf{x}$ with respect to time. Let's do this first. We get:

$$
\mathbf{x}=\mathbf{v}e^{i\omega t}\\[1em]
\mathbf{\dot{x}}=i\omega\mathbf{v}e^{i\omega t}\\[1em]
\mathbf{\ddot{x}}=i^2\omega^2\mathbf{v}e^{i\omega t}
$$

Since we are using complex numbers here, we know that $i=\sqrt{-1}$ by definition, and so we have $i^2=-1$ by definition as well. Using this, we can simply our second-order derivative $\mathbf{\ddot{x}}$ to:

$$
\mathbf{\ddot{x}}=-\omega^2\mathbf{v}e^{i\omega t}
$$

If we insert this now into Eq.(\ref{eq:rigid-body-motion-general}), and assuming that we have $f=0$, then we get:

$$
\mathbf{M}\left(-\omega^2\mathbf{v}e^{i\omega t}\right)+\mathbf{K}\left(\mathbf{v}e^{i\omega t}\right)=0
$$

We can factor out the term $\exp({i\omega t})$ which gives us:

$$
\left(-\omega^2\mathbf{M}\mathbf{v}+\mathbf{K}\mathbf{v}\right)e^{i\omega t}=0
$$

Since the right-hand side is zero, we can divide by $\exp{i\omega t}$ to get rid of it (and, with that, we loose any trace of complex numbers, hooray!), which gives us:

$$
-\omega^2\mathbf{M}\mathbf{v}+\mathbf{K}\mathbf{v}=0
$$

We rewrite this equation as:

$$
\mathbf{K}\mathbf{v} = \omega^2\mathbf{M}\mathbf{v}
$$

We now define $\lambda=\omega^2$ and obtain:

$$
\mathbf{K}\mathbf{v} = \lambda\mathbf{M}\mathbf{v}
$$

Since $\omega$ is the angular frequency, it is probably of great importance and so we want to know the values of $\lambda=\omega^2$. Here, $\lambda$ can have many values and these are our eigenvalues. For example, the smallest eigenvalue may tell us something about the dominant frequencies that will be relevant for studying the stability of ships.

So, how do we compute them? There are different ways of doing it, and they are not really of great importance here, but one approach could be to multiply by the inverse of the $\mathbf{M}$ matrix which would give us:

$$
\mathbf{M}^{-1}\mathbf{K}\mathbf{v} = \lambda\mathbf{v}
$$

Let's look at the left-hand side in some more detail. We have the multiplication of $\mathbf{Kv}$. Here, $\mathbf{K}$ is a matrix, while $\mathbf{v}$ is the mode vector. Thus, the product $\mathbf{Kv}$ is a vector. Let's call this vector $\mathbf{b}$. Thus, we now have to evaluate the product of $\mathbf{M}^{-1}\mathbf{b}$. Let's say the solution of this is given by the vector $\mathbf{x}$. If we replace the matrix $\mathbf{M}$ with $\mathbf{A}$ (we just change the symbol, not the meanting of the matrix itself), then we are solving:

$$
\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}
$$

This is the equation we have to solve if we have the general linear system of equations of the form $\mathbf{Ax}=\mathbf{b}$, which we have now seen a few times already. So, you see, Krylov was interested in ship dynamics, but the fundamental problem he had to solve was $\mathbf{Ax}=\mathbf{b}$, or rather, how to form the product of $\mathbf{A}^{-1}\mathbf{b}$.

Here is the trick, or what I have referred to as *cheating*. Krylov didn't form the inverse of $\mathbf{A}^{-1}$ (or $\mathbf{M}^{-1}$). No, instead, he formed a subspace of matrix vector products of the form

$$
\mathcal{K}_r(\mathbf{A}, \mathbf{v})=\text{span}\{\mathbf{v}, \mathbf{A}\mathbf{v}, \mathbf{A}^2\mathbf{v}, \mathbf{A}^3\mathbf{v}, ..., \mathbf{A}^{r-1}\mathbf{v}\}
$$

All of this multiplications where much easier to perform, but, it turns out that these matrices he was investigating behaved roughly the same way as the original matrix $\mathbf{M}$ (i.e. he was making burgers consisting of only a bun and a patty, instead of burgers with additional flavours like sauce, onions, tomato, etc.)

This was his trick, and many methods in the field of linear algebra were developed ontop of his subspace idea. For example, you may come across methods like the Lanczos or Arnolid method; these methods are build ontop of the Krylov subspace and allow us to approximate with great ease and acceptable accuracy the dominant eigenvalues in a matrix.

Why approximate? Well, the algorithm's themselves are exact, so, if we applied them to the original matrix $\mathbf{M}$, for example, we would get all eigenvalues computed exactly. But, that would be computationally inefficient (and, infact, in most cases also not possible). So, instead, these algorithms make use of the Krylov subspace, which make problems much smaller and therefore easier to compute.

You can think of Lanczos and Arnoldi as two burger preparing experts at Burger King; if they only use a bun and a patty to make any burger on the menu, well, then they are going to be extremely fast at their work. Regardless of the burger that is ordered, their preparation time is always the same!

If I was the acting manager, I would take notice at how quickly Lanczos and Arnoldi are working and, heck, I may even promote them to supervisors as they are clearly doing something right! Well, but if I were the acting manager and then eat one of the burgers they have prepared, I would realise that yes, they are fast, but not accurate; they have approximated a burger, but it is not exactly what I would have expected!

Lanczos and Arnoldi give us approximation of the eigenvalues based on the simplified Krylov subspace which still captures the essence of the original problem. They are fast, but approximate. For our applications, that is fine, we don't care a great deal about all eigenvalues; typically, we only care about the dominat eigenvalues, that is, the smallest and largest eigenvalues.

Without looking at the algorithms in detail here, Arnoldi is a general algorithm we can always apply to a matrix to find its eigenvalues. Lanczos is a specialisation that works only for symmetric matrices.

So, in summary, the Krylov subspace makes the computation of certain quantities of interest, like the Eigenvalues, really easy. On top of that, a linear system of equation solver for $\mathbf{Ax}=\mathbf{b}$ can be constructed by exploiting the Krylov subspace idea. Doing so will lead to a number of algorithms that are simply known as Krylov subspace methods. We'll get to those later, but for now, let's look at why eigenvalues are important to us!

### The role of eigenvalues in linear algebra

First of all, eigenvalues are not really of importance to us in the sense that we want to compute them, but more in the sense that knowing the eigenvalues of certain matrices will typically tell us something about the convergence speed.

The first important property is the so-called spectral radius. The spectral radius of a matrix $\mathbf{A}$ is the largest eigenvalue of $\mathbf{A}$. We use the letter $\rho$ to define the spectral radius and write its definition as:

$$
\rho(\mathbf{A}) = \text{max}\{|\lambda_i|\}
\tag{eq:spectral-radius}
$$

Here, we may have many eigenvalues, but we only care about the largest eigenvalue. So, for example, if we had:

$$
\lambda=(3,-2,1,5,-1)
$$

Then $\rho(\mathbf{A})=5$. Later, when we deal with coefficient matrices $\mathbf{A}$, we will be interested in the spectral radius of these coeffcient matrices. As it turns out, we have the following conditions:

$$
\rho(\mathbf{A})=
\begin{cases}
\lt 1\qquad &\text{Convergence}\\
1\qquad &\text{Stagnation}\\
\gt 1\qquad &\text{Divergence}
\end{cases}
$$

Thus, if the largest eigenvalue is less than 1, that is, our sepctral radius is $\rho(\mathbf{A})\lt 1$, our iterative system will converge. It is greater than 1, it will diverge, and, if it is 1, then we have neither convergence nor divergence.

Furthermore, the smaller the largest eigenvalue is, the faster the convergence will become. So, ideally, we want to construct coefficient matrices that have extremely small eigenvalues.

Another property of our coefficient matrix is linked to the so-called condition number. As long as our matrix is symmetric positive definite (a property we will look at in the next section), we can compute the condition number of our coefficient matrix as the ratio if its largest and smallest eigenvalues. This is given as:

$$
\kappa(\mathbf{A})=\frac{\lambda_{max}}{\lambda_{min}}
\tag{eq:condition-number}
$$

Why is this condition number important? It tells us something about convergence again. The smaller the condition number is (the closer the minimum and maximum eigenvalues are together), the faster the convergence becomes. Thus, we would like to construct coefficient matrices which have very narrowly spread eigenvalues to have condition numbers close to 1. This will give us the fastest convergence.

This will then later lead to the idea of preconditioners, where we replace our original coefficient matrix by an equivalent matrix that will give us the same results, but one that has a much smaller condition number.

So, in summary, we don't necessarily need to compute eigenvalues when we solve a linear system of equations of the form $\mathbf{Ax}=\mathbf{b}$, but we need them if we want to analyse the convergence rate of our algorithms to approximate the solution $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$. Thankfully, all of that has been already done, so we have a good idea of which method works well and which convergeces the fastest, so we don't have to do this everytime we solve a linear system of equations!

### Definite matrices

All matrices we will generate as part of our discretisation process will necessarily be square matrices. The number of rows and columns is the same, and it is set to the number of cells, or vertices, depending on where we store our solution variables.

When we talk about the definiteness in a linear algebra sense, we look at our matrix $\mathbf{A}$, and what it would do to a vector $\mathbf{x}$ if we multiplied it in the following form:

$$
\mathbf{x}^T\mathbf{Ax}
$$

The product of $\mathbf{Ax}$ gives a column vector. If we multiply that by a row vector, i.e. $\mathbf{x}^T$, then we get a scalar. Based on the values we obtain, we can characterise our matrix into the following forms:

A matrix is said to be *positive definite* if:

$$
\mathbf{x}^T\mathbf{Ax}\gt 0,\qquad \mathbf{x}\ne\mathbf{0}
$$

That is, any vector multiplied in the above form with a vector containing not only zeros, will produce a scalar greater than zero. If a matrix is positive definite, than all of its eigenvalues are real and greater than zero. This is an important property we will later use.

A matrix is said to be *negative definite* if:

$$
\mathbf{x}^T\mathbf{Ax}\lt 0,\qquad \mathbf{x}\ne\mathbf{0}
$$

A matrix is said to be *positive semidefinite* if:

$$
\mathbf{x}^T\mathbf{Ax}\ge 0
$$

And, finally, a matrix is said to be *negative semidefinite* if:

$$
\mathbf{x}^T\mathbf{Ax}\le 0
$$

Some matrices may be symmetric, that is, taking the transpose of the matrix does not change the matrix itself and we have $\mathbf{A}^T=\mathbf{A}$. When we talk about definite matrices, we typically assume the matrix $\mathbf{A}$ to be symmetric. This isn't the case for all coefficient matrices that we will obtain; the upwind scheme, for example, will destroy symmetry.

But, equations without convection (for example, the pressure Poisson equation) can be created in a symmetric form. It depends on how we impose boundary conditions, but it is possible to obtain a positive definite matrix during the discretisation.

This is important, as some algorithms assume Symmetric, Positive Definite (or SPD) matrices when we deal with linear algebra solvers that try to find an approximation to $\mathbf{Ax}=\mathbf{b}$. If these matrices are not SPD, then the algorithms will break down. For example, the Conjugate Gradient (CG) method requires SPD matrices, and it will stop working for matrices that look differently.

I don't know how many times I have implemented a matrix that looked alright (SPD), which I then tried to solve with a CG algorithm, only to then realise that the boundary conditions, or something else, were breaking the symmetry. For these cases, we have generalisation of the CG algorithm, like the Bi-Conjugate Gradient Stabilised (BiCGStab) method, which can work with non-symmetric matrices, but I am getting ahead of myself.

The important take away here is that some matrices are positive definite (and symmetric), and if they are, we can write bespoke algorithms for them for fast convergence. And with that, I'd say we have looked at all the matrix properties that will be of importance to us. Let us now turn to solution algorithms to solve matrices.

## Direct methods

A direct method is one where we seek to solve $\mathbf{Ax}=\mathbf{b}$ exactly. We can either create $\mathbf{A}^{-1}$ directly, or manipulate $\mathbf{A}$ in such a way that computing the unknown vector becomes trivial.

In our introductory example, we obtained the following coefficient matrix for our explicit time integration:

$$
\mathbf{A}=
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
0 & 1/\Delta t & 0 & 0 & 0 \\[1em]
0 & 0 & 1/\Delta t & 0 & 0 \\[1em]
0 & 0 & 0 & 1/\Delta t & 0 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
$$

Since we only have diagonal entries, we can easily invert this and write:

$$
\mathbf{A}^{-1}=
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\[1em]
0 & \Delta t & 0 & 0 & 0 \\[1em]
0 & 0 & \Delta t & 0 & 0 \\[1em]
0 & 0 & 0 & \Delta t & 0 \\[1em]
0 & 0 & 0 & 0 & 1 \\[1em]
\end{bmatrix}
$$

Now we can calculate $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$ with ease. In fact, whenever we look at an explicit discretisation, we have an equation of the form:

$$
x^{n+1}_i=x^n_i + \Delta t(...)
$$

Here, $(...)$ denotes any term that were in the discretised equation that we brought onto the right-hand side. But we can see that we multiply everything on the right-hand side here with $\Delta t$; this is $\mathbf{A}^{-1}$.

For implicit time integrations, $\mathbf{A}^{-1}$ isn't easily invertable, as we already saw, and so now we have to start looking at solving $\mathbf{Ax}=\mathbf{b}$ somehow. In this section, we will look at direct methods, though these are usually not used in production codes as they are too slow. Some are better than others, and some of these we will come across later again when we deal with preconditioning, so even if we don't use direct methods in CFD directly, we use them in modified forms.

### Gauss Elimination

The first, and probably most well know algorithm is the Gaussian elimination process. You'll likely have come across that already in high-school, and perhaps you have already forgotten it. If you have done CFD for some time, you may say that you have never used the Gaussian elimination in any CFD solver. That is true, there are good reasons, as we will see, that we don't use the Gaussian elimination in practice, but, this algorithm forms the basis for a lot of other algorithms that we do use.

Thus, let us have a look how Gaussian elimination works, so that we will be comfortable when we need it in subsequent sections. Our goal is to solve $\mathbf{Ax}=\mathbf{b}$, but without forming the inverse of $\mathbf{A}$. If we could easily and cheaply computed $\mathbf{A}^{-1}$, then we would not need this article, i.e. everything we discuss in this article discusses how to solve $\mathbf{Ax}=\mathbf{b}$ without $\mathbf{A}^{-1}$.

The way the Gaussian elimination achieves that is by manipulating $\mathbf{A}$ and $\mathbf{b}$ until it becomes trivial to solve the system. The goal here is to make $\mathbf{A}$ an upper triangular matrix, so, if we have $\mathbf{A}$ given as:

$$
\mathbf{A}=
\begin{bmatrix}
a & b & c \\
d & e & f \\
g & h & i
\end{bmatrix}
$$

Then, we want to bring this into a modified (upper triangular matrix) form as:

$$
\mathbf{A}=
\begin{bmatrix}
a & b & c \\
0 & e' & f' \\
0 & 0 & i'
\end{bmatrix}
$$

Here, we have not changed row 1, but all rows below will now have modified coefficients, as indicated by $e'$, $f'$, and $i'$, respectively. If the matrix is given in this form, we can write the full system of $\mathbf{Ax}=\mathbf{b}$ as:

$$
\begin{bmatrix}
a & b & c \\
0 & e' & f' \\
0 & 0 & i'
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}=
\begin{bmatrix}
b_1 \\
b'_2 \\
b'_3
\end{bmatrix}
\tag{eq:gaussian-elimination-idea-example}
$$

Notice here that if we have to change row 2 and 3 in the $\mathbf{A}$ matrix, we also have to make changes in the $\mathbf{b}$ matrix. If we have done that, and we look at the last row in Eq.(\ref{eq:gaussian-elimination-idea-example}), we can write out the equation for the last row as:

$$
i'\cdot x_3=b'_3
$$

Solving this for $x_3$ and we get:

$$
x_3 = \frac{b'_3}{i'}
$$

If we now look at the second row in Eq.(\ref{eq:gaussian-elimination-idea-example}), we can write this as:

$$
e'\cdot x_2 + f'\cdot x_3 = b'_2
$$

Since we already have computed $x_3$ in the previous step, we can write this equation as:

$$
e'\cdot x_2 + f'\cdot \frac{b'_3}{i'} = b'_2
$$

Now, this equation only contains a single unknown, and so we can solve for $x_2$ as:

$$
x_2 = \frac{1}{e'}\left(b'_2 - f'\cdot \frac{b'_3}{i'}\right)
$$

We can do the same now for the first row to find $x_1$, and in this way, we have a solution for the unknown vector $\mathbf{x}$ from our linear system of equations $\mathbf{Ax}=\mathbf{b}$ without having to form $\mathbf{A}^{-1}$. Since we obtain values in the $\mathbf{x}$ vector from the back to the front (starting at $x_3$ in our example and going to the first entry), we call this step the backward substitution.

OK, so how do we get both $\mathbf{A}$ and $\mathbf{b}$ into this modified form? Well, this is best illustrated with an example. Let's say we want to solve the following system:

$$
\begin{bmatrix}
-2 & 1 & 0 \\
1 & -2 & 1 \\
0 & 1 & -2
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}=
\begin{bmatrix}
-3 \\
3 \\
-5
\end{bmatrix}
$$

Our matrix $\mathbf{A}$ here could have been obtained, for example, for the discretisation of the pressure Poisson equation, i.e. it has the same entries on the diagonal and off-diagonals. It is a sparse matrix, but with only 3 points, this sparsity does not really show. In any case, let's say we want to solve this equation now using Gaussian elimination.

We want to bring $\mathbf{A}$ into an upper triangular form. First of all, we write $\mathbf{A}$ and $\mathbf{b}$ into a compact form as:

$$
\left[\begin{array}{c c c | c}
-2 & 1 & 0 & -3 \\
1 & -2 & 1 & 3 \\
0 & 1 & -2 & -5
\end{array}\right]
$$

We leave the first row as it is and go to the second row. If we want to get an upper triangular matrix, that is, one where all the entries below the diagonal of $\mathbf{A}$ are zero, we need to eliminate all of the values below the diagonal. We achieve that by adding, or subtracting, multiples of the first row.

So, if we want to produce zeros in the first column, we have to subtract multiples of the first row from all rows below. If we want to produce zeros in the second column, then we have to subtract multiples of the second row from all rows belot it, and so on.

So, let's look at the second row. We want to get rid of the first entry (first column) in the second row, which is 1. I said that if we want to get rid of elements in the first column, we have to subtract multiples of the first row. We could write this in pseudo maths as:

$$
\text{Row 2} - \text{multiplier}\cdot\text{Row 1}
$$

I like to call this also business math, or finance math. If you have ever seen text books on finance who claim "we can do math as well", you know what I mean. In any case. We want to get rid of the entry in the first column of the second row, and so we can write the equation

$$
1 - m\cdot (-2) = 0
$$

Here, 1 is the first entry from row 2, and we subtract the first entry from row 1, multiplied by some value so that this becomes zero. So, if we take $m=-0.5$, we get:

$$
1 - (-0.5)\cdot (-2) = 1 - 1 = 0
$$

So, we have to now multiply all entries in the first row by $m=-0.5$ and subtract that from the second row. In our compact system, we can write this as:

$$
\left[\begin{array}{c c c | c}
-2 & 1 & 0 & -3 \\
0 & -1.5 & 1 & 1.5 \\
0 & 1 & -2 & -5
\end{array}\right]
$$

We have to do the same for all rows below it until we have zeros everywhere, expect for row 1. As it turns out, the last row already contains a zero here, so we don't have to do anything, and we can move onto the next row.

Now, we are in row 2, and so we want to produce zeros in the second column below row 2. This means, the second entry in the third row, which has a value of 1, needs to be eliminated. So, the question again becomes, what do I have to multiply my second row with so that when I subtract it from the third row, the entry in the second column will become zero for the third row.

We can write this again as an equation:

$$
1 - m\cdot (-1.5) = 0
$$

If we pick $m=-2/3$, then we satisfy the above equation. To verify this, let's write this out:

$$
1 - \left(-\frac{2}{3}\right)\cdot (-1.5) = 1 - \left(-\frac{2}{3}\right)\cdot \left(-\frac{3}{2}\right) = 1 - \frac{6}{6} = 1 - 1 = 0
$$

OK, so let's multiply the second row by $m-2/3$ and subtract that from the third row. This will result in:

$$
\left[\begin{array}{c c c | c}
-2 & 1 & 0 & -3 \\
0 & -1.5 & 1 & 1.5 \\
0 & 0 & -4/3 & -4
\end{array}\right]
$$

As we can see, by eliminating entries in the first and second column below the diagonal, we have obtained an upper triangual matrix. So let's write our modified system $\mathbf{Ax}=\mathbf{b}$ in full:

$$
\begin{bmatrix}
-2 & 1 & 0 \\
0 & -1.5 & 1 \\
0 & 0 & -4/3
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}=
\begin{bmatrix}
-3 \\
1.5 \\
-4
\end{bmatrix}
$$

Using the backward substitution, we can now find the values of $\mathbf{x}$. Let's do that. The last row can be written as:

$$
-\frac{4}{3}x_3 = -4
$$

Solving this for $x_3$ gives us:

$$
x_3 = -4\cdot \left(-\frac{3}{4}\right) = \frac{12}{4}= 3
$$

So, we know the first entry in $\mathbf{x}$. We can use that knowledge and write the equation in the second row as:

$$
-1.5\cdot x_2 + x_3 = 1.5
$$

Inserting $x_3 = 3$ and solving for $x_2$, we get:

$$
-1.5\cdot x_2 + 3 = 1.5
-1.5\cdot x_2 = 1.5 - 3 = -1.5
x_2 = \frac{-1.5}{-1.5} = 1
$$

And so, we know that $x_2=1$. Now that we know this value, we can go to the first row and write it out as:

$$
-2\cdot x_1 + x_2 = -3
$$

Inserting $x_2 = 1$, we have:

$$
-2\cdot x_1 + 1 = -3
-2\cdot x_1 = -3 -1 = -4
x_1 = \frac{-4}{-2}=2
$$

Thus, with $x_1=2$ determined, we have found the solution to $\mathbf{x}$ as $\mathbf{x}=(2,1,3)^T$. We can verify that by computing the product $\mathbf{Ax}$:

$$
\begin{bmatrix}
-2 & 1 & 0 \\
1 & -2 & 1 \\
0 & 1 & -2
\end{bmatrix}
\begin{bmatrix}
2 \\
1 \\
3
\end{bmatrix}=
\begin{bmatrix}
-3 \\
3 \\
-5
\end{bmatrix}
$$

Comparing the result with our right-hand side vector $\mathbf{b}$, given as $\mathbf{b}=[-3, 3 -5]^T$, we can see that the solution we found for $\mathbf{x}$ produces the same right-hand side vector, and all of that without ever having to compute $\mathbf{A}^{-1}$.

And that is the Gaussian elimination. We get an exact solution for all entries in $\mathbf{x}$, and so, you may say, great, this method is solid and we can use that for every linear system of equation of the form $\mathbf{Ax}=\mathbf{b}$. Well, yes, we can, but it turns out, the Gaussian elimination is just too expensive to perform.

Let's look at the cost. Let's generalise our matrix and say it has $n$ entries. In the example we looked at, we had $n=3$. First, we have to loop through all rows (except the first one), which will take us $n$ operations. Then, for each row, we have to go through each column and subtract the first row (if we want to write zeros into the first column, or the N-th row, if we want to write zeros into the N-th column).

So, we loop over each row, which costs $n$ operation, but for each row, we have to loop over $n$ columns, so the total computational cost for that is of the order of $n^2$. This is only to bring the matrix into upper triangular form, but then, we also have to perform the backward substitution. This requires us to loop over $n$ rows again, increasing our computational cost to $n^3$.

This is a problem. Why? Well, let's add some numbers to make this more real. Let's say we have $n=100$ and we implement the Gaussian elimination algorithm. We solve our system, and for the sake of argument, let's say we time how long it takes for our algorithm to perform this computation and we get 1 second as the answer.

So, if we use $n=100$, and we have a cost of $n^3=100^3=100\cdot 100\cdot 100 = 1,000,000$, then we can say, that $1,000,000$ operations are equivalent to 1 second of computational cost.

OK, so instead of having a matrix with $n=100$ entries, let's say we have $n=1000$. The computational cost would be $n^3=1,000\cdot 1,000\cdot 1,000=1,000,000,000$. If we said that $n=100$ required 1 million operations, which took 1 second to solve, then we can say that $n=1,000$ will require 1 billion operations. If we divide 1 billion by 1 million, we get 1000 as a result, and so we can say that by increasing the number of points by a factor of 10, we increased the computational cost by a factor of 1000.

So, insteado of having to wait for 1 second, we have to wait for 1,000 seconds now. And now imagine what happens if we start to use more realistic numbers. In reality, we use millions, even hundreds of millions of elements in our matrix. The computational cost to solve these linear system, i.e. $\mathbf{Ax}=\mathbf{b}$, with the Gaussian elimination is so high that we don't even bother.

So, Gaussian elimination is a nice example that we can go through on paper, or in a classroom, and we get a solution for $\mathbf{x}$, but in a real CFD code, we would never use it due to its prohibitively high computational cost. Having said that, we can generalise this procedure somewhat, which results in the so-called LU decomposition, or LU factorisation.

This decomposition is based entirely on the Gaussian elimination, and we do use this in real CFD codes, due to some properties we can exploit in the LU decomposition. So, even if we don't use Gaussian elimination directly in CFD codes, it still shows up, in disguise, in other methods. So, let's have a look then at what this LU decomposition is.

### Lower-Upper (LU) Decomposition

Indeed, the LU decomposition is just the Gaussian elimination with some structure, and to see why it works, we have to look at some matrix properties first.

If I have a matrix $\mathbf{A}$, then I can decompose it into separate matrices and add them together if I want. For example, if we take the same matrix $\mathbf{A}$ that we can had before:

$$
\mathbf{A}=
\begin{bmatrix}
-2 & 1 & 0 \\
1 & -2 & 1 \\
0 & 1 & -2
\end{bmatrix}
$$

Then, we can decompose this matrix into several other matrices as:

$$
\mathbf{A}=\mathbf{B} + \mathbf{C} + \mathbf{D}=
\begin{bmatrix}
-2 & 0 & 0 \\
1 & 0 & 0 \\
0 & 0 & 0
\end{bmatrix}+
\begin{bmatrix}
0 & 1 & 0 \\
0 & -2 & 0 \\
0 & 1 & 0
\end{bmatrix}+
\begin{bmatrix}
0 & 0 & 0 \\
0 & 0 & 1 \\
0 & 0 & -2
\end{bmatrix}
$$

The idea behind the LU decomposition, or LU factorisation, is that we split the matrix into a lower and upper triangular matrix. So, we may, naively assume that the LU decomposition is just the decomposition of $\mathbf{A}$ into its lower and upper triangular matrix as:

$$
\mathbf{A}=\mathbf{L}+\mathbf{U}=
\begin{bmatrix}
-2 & 1 & 0 \\
0 & -2 & 1 \\
0 & 0 & -2
\end{bmatrix}+
\begin{bmatrix}
0 & 0 & 0 \\
1 & 0 & 0 \\
0 & 1 & 0
\end{bmatrix}
$$

While this decomposition is correct, it is not what we mean by the LU decomposition. Instead of decomposing the matrix $\mathbf{A}$ into its lower and upper triangular matrix which can be *added* back together to result in the original matrix $\mathbf{A}$, we want to find matrices that can be *multiplied* together to give $\mathbf{A}$. So, the LU decomposition tries to decompose the matrix $\mathbf{A}$ as:

$$
\mathbf{A}=\mathbf{LU}
$$

The starting point for the LU decomposition is to write the product of the matrix $\mathbf{A}$ with the identity matrix $\mathbf{I}$ as:

$$
\mathbf{A}=\mathbf{IA}
\begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
-2 & 1 & 0 \\
1 & -2 & 1 \\
0 & 1 & -2
\end{bmatrix}
$$

We want to manipulate both matrices here so that multiplied together they still result in $\mathbf{A}$, but at the same time, both matrices form a lower (L) and upper (U) triangular matrix. We saw in the previous section that when we do Gaussian elimination, that we already get the upper triangular matrix. So we know how to form this matrix, and this is what we will be doing in a second. But what about the lower triangular matrix L?

Well, this matrix will simply store the coefficients we use to multiply the equations with in order to zero out a column. So, we keep track of the multiplication factors we use in the Gaussian elimination. Let's do that with our matrix $\mathbf{A}$ and see how this works in practice.

When we started the Gaussian elimination process, we said that we needed to multiply the first row by $m-0.5$ and then subtract this first row from the second. So, we will do that now again, and record the result in the second row in what will become our upper triangular matrix (what currently is $\mathbf{A}$ above) and we will also store $m-0.5$ in what currently is $mathbf{I}$ (and what will become our lower triangular matrix $\mathbf{L}$)

We store $m-0.5$ in the second row and first column, i.e. the location where we want to zero out the matrix $\mathbf{A}$. If we apply that for the second row, then we get:

$$
\mathbf{A}=
\begin{bmatrix}
1 & 0 & 0 \\
-0.5 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
-2 & 1 & 0 \\
0 & -1.5 & 1 \\
0 & 1 & -2
\end{bmatrix}
$$

The third row still has a zero in the first column, so we do not need to change that. Next, we go to the third row, to the second column, and we want to get rid of the entry here as well, so that we get an upper triangular structure for the matrix on the right above. We said that we can subtract the second row from the third, if we multiply the second row by $m=-2/3$. Then we can record the value of $m$ in the lower triangular matrix again and modify our upper triangular matrix and get:

$$
\mathbf{A}=
\begin{bmatrix}
1 & 0 & 0 \\
-0.5 & 1 & 0 \\
0 & -2/3 & 1
\end{bmatrix}
\begin{bmatrix}
-2 & 1 & 0 \\
0 & -1.5 & 1 \\
0 & 0 & -4/3
\end{bmatrix}
$$

And that's it! The first matrix is our lower triangular matrix $\mathbf{L}$, while the second matrix is the upper triangular matrix $\mathbf{U}$. Again, this upper triangular matrix is obtained from the Gaussian elimination, while we store the coefficients of multiplications of the different rows in the lower triangular matrix $\mathbf{L}$. To see that this decomposition is valid, let's multiply $\mathbf{L}$ and $\mathbf{U}$ together to verify that this results again in $\mathbf{A}$.

For example, to obtain the first row of $\mathbf{A}$, we need to multiply the elements of the first row of $\mathbf{L}$ with the columns of $\mathbf{U}$ and then add the results together. For the first row and first column, we get:

$$
a_{11}= 1\cdot (-2) + 0 \cdot 0 + 0\cdot 0 = -2
$$

For the second entry in the first row, we multiply the first row of $\mathbf{L}$ with the second column of $\mathbf{U}$ and get:

$$
a_{12} = 1 \cdot 1 + 0\cdot (-1.5) + 0\cdot 0 = 1
$$

And, for the third entry, we multiply the first row of $\mathbf{L}$ with the third column of $\mathbf{U}$ and get

$$
a_{13} = 1\cdot 0 + 0\cdot 1 + 0\cdot (-4/3) = 0
$$

In the same way, we can find the second and third row in $\mathbf{A}$ and get:

$$
\mathbf{LU}=
\begin{bmatrix}
1 & 0 & 0 \\
-0.5 & 1 & 0 \\
0 & -2/3 & 1
\end{bmatrix}
\begin{bmatrix}
-2 & 1 & 0 \\
0 & -1.5 & 1 \\
0 & 0 & -4/3
\end{bmatrix}=
\begin{bmatrix}
-2 & 1 & 0 \\
1 & -2 & 1 \\
0 & 1 & -2
\end{bmatrix}
$$

Indeed, this type of multiplication results again in the original coefficient matrix $\mathbf{A}$. So, you may then rightfully ask, what's the point? We said before that Gaussian elimination isn't great from a computational cost point of view, so why should we bother with the LU decomposition then, which is just the Gaussian elimination, but written in a slightly different form?

Well, the main advantage of the LU decomposition is that we can reuse it, while the Gaussian elimination cannot be reused. Did you notice that we used the right-hand side vector $\mathbf{b}$ in the Gaussian elimination, which was completely absent in the LU decomposition? This means that if the right-hand side vector $\mathbf{b}$ changes, we have to do the Gaussian elimination process again, whereas the LU decomposition can be reused even if the right-hand side vector changes.

This is one of its main advantages. To see this, let's look at how we would compute the unknown solution vector $\mathbf{x}$ in a linear system of equation $\mathbf{Ax}=\mathbf{b}$ with the LU decomposition. In the Gaussian elimination, we used a backward substitution. Using the LU decomposition, we have to use a forward and backward substitution. Formally, we can write this as:

$$
\mathbf{Ly}=\mathbf{b}
$$

This is the forward substitution. The backward substitution can be written as:

$$
\mathbf{Ux}=\mathbf{y}
$$

The vector $\mathbf{y}$ is an intermediate result that helps us to solve for the vector we are actually interested in, i.e. $\mathbf{x}$. In our example in the Gaussian eliminations ection, we used the right-hand side vector $\mathbf{b}=(-3,3,-5)^T$, so, let's use that again, as well as our lower and upper triangular matrices $\mathbf{L}$ and $\mathbf{U}$ to find the solution vector $\mathbf{x}$.

First, we have to do our forward substitution. This is:

$$
\mathbf{Ly}=\mathbf{b}\\[1em]
\begin{bmatrix}
1 & 0 & 0 \\
-0.5 & 1 & 0 \\
0 & -2/3 & 1
\end{bmatrix}
\begin{bmatrix}
y_1 \\
y_2 \\
y_3
\end{bmatrix}=
\begin{bmatrix}
-3 \\
3 \\
-5
\end{bmatrix}
$$

It is called a forward substitution because we solve for $\mathbf{y}$ by obtaining each value within this vector by going from the front to the back. For example, $y_1$ can be obtained as:

$$
1\cdot y_1 = -3 \\[1em]
y_1 = -3
$$

We can obtain $y_2$ from the second row as:

$$
-0.5\cdot y_1 + 1\cdot y_2 = 3 \\[1em]
-0.5\cdot (-3) + 1\cdot y_2 = 3 \\[1em]
1.5 + y_2 = 3 \\[1em]
y_2 = 1.5 \\[1em]
$$

And, finally, from row 3, we can obtain $y_3$ as:

$$
-2/3\cdot y_2 + 1\cdot y_3 = -5\\[1em]
-2/3\cdot 1.5 + 1\cdot y_3 = -5\\[1em]
-1 + y_3 = -5\\[1em]
y_3 = -4\\[1em]
$$

Thus, we have found the intermediate result $\mathbf{y}=(-3, 1.5, -4)^T$. This was the modified right-hand side vector when we constructed the Gaussian elimination. With this vector $\mathbf{y}$, we can now solve the backward substitution as:

$$
\mathbf{Ux}=\mathbf{y}\\[1em]
\begin{bmatrix}
-2 & 1 & 0 \\
0 & -1.5 & 1 \\
0 & 0 & -4/3
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3
\end{bmatrix}=
\begin{bmatrix}
-3 \\
1.5 \\
-4
\end{bmatrix}
$$

Well, we can do that, or recognised that we have already done that in the Gaussian elimination. As I have mentioned, the right-hand side vector $\mathbf{y}$ is the same as the modified right-hand side vector $\mathbf{b}'$ from the Gaussian elimination process we saw before, and the upper triangular matrix $\mathbf{U}$ is also what we obtained in the Gaussian elimination process.

Once we had this matrix, we were then able to compute the values for $\mathbf{x}$ from the back, hence this is the backward substitution process. For completeness, we found the solution vector to be $\mathbf{x}=(2,1,3)^T$. If you want to verify that, we have to write out the equation for the last row, and we get:

$$
-4/3\cdot x_3 = -4 \\[1em]
x_3 = -4\cdot (-3/4) \\[1em]
x_3 = 3
$$

We can do this for the second and thrid row as well to confirm that $x_2=1$ and $x_1=2$. OK, so we states thus far that the LU decomposition is just the Gaussian elimination in disguise, we store the Gaussian elimination in a slightly different form, which allows us to solve for the unknown solution vector $\mathbf{x}$ through a forward and backward substitution. What gives?

I said that once $\mathbf{L}$ and $\mathbf{U}$ are found and $\mathbf{A}$ doesn't change, then we can reuse these matrices to compute the unknown solution vector $\mathbf{x}$ in $\mathbf{Ax}=\mathbf{b}$, that is, the right-hand side vector $\mathbf{b}$ can change, but $\mathbf{A}$ is not. This is the case for a few examples in CFD, for example, the pressure Poisson equation typically has a constant coefficient matrix, and here an LU decomposition would make sense.

We still have the cost of having to perform the Gaussian elimination once, and for large systems, this becomes quite expensive. So, instead, when we deal with the LU decomposition in CFD, we use the so-called incomplete LU decomposition, or ILU in short. We will see this later when we talk about preconditioners, where the LU decomposition is typically used.

As it turns out, performing an approximate LU decomposition can drastically reduce the computational cost of the original Gaussian elimination, so much so that it becomes feasible again for CFD calculations. However, it is approximate, and by far not accurate, so we have to couple it with another method to find an accurate solution. Again, we will talk about that in more detail when we deal with preconditioners.

What about the cost for the LU decomposition? Well, we still have the cost of $n^3$ to find both $\mathbf{L}$ and $\mathbf{U}$. But, once these are found, and assuming $\mathbf{A}$ isn't changing, we can reuse these matrices and only perform the forward and backward substitution. So, how expensive are these?

In the forward substitution, we have to go through each row of the $\mathbf{L}$ matrix, which costs $n$ operations. We then have to multiply each row with values in $\mathbf{y}$, i.e. in each column, which has a cost of $n$ as well (we have to loop over each column for each row). Thus, the cost for the forward substitution is $n^2$.

The same arguments can be made for the backward substitution, which also has cost $n^2$. Performing two operations that cost $n^2$ is still better than performing one operation with a cost of $n^3$, so the LU decomposition wins here over the Gaussian elimination. For the incomplete LU decomposition, we will see that the cost becomes more like $n$ for the forward and backward substitution, so we can even further reduce the cost, which makes it very lucrative.

So, we saw how the LU decomposition is just the Gaussian elimination in disguise, and while we will not use it in the form we have discussed in this section, we will later form an approximation, that will allow us to reduce computational costs substantially. Textbooks on CFD will have you believe that we don't use direct methods, or the Gauss elimination, but we do, just not in the textbook form.

I want to look at one more powerful direct methods that does actually have a place in CFD and, for the right application, may be the best direct solver you can implement. This will be discussed in the next section.

### Thomas and his little algorithm

Llewellyn Hilleth Thomas was an interesting man. Born in 1903 in London, 57 day before the first powered flight by the wright brothers in the USA, he was born right around the time aviation was born. I would like to think that he became known as *baby jesus of aviation*, but Llewellyn had different plans (and yes, I have no idea how to pronounce that name either, even after 15 years of self-imposed exile in the UK, I can still not pronounce welsh names. But, brits can't either, so I am not too bothered ...). Let us call him Thomas from now on.

Thomas went on to have a successful career in Physics. While he turned his back to the field of aviation, he did experience the birth of quantum physics in his twenties, and back then, with the absence of Tik Tok, Taylor Swift, and Snapchat, well, quantum physics must have been *the most* sexy thing around, and so he threw himself at it.

He developed a relativistic correction for the spin of particles, known as the [Thomas precision](https://en.wikipedia.org/wiki/Thomas_precession), as well as a theory for electronic structure of many-body systems, known as the [Thomas-Fermi (TF) model](https://en.wikipedia.org/wiki/Thomas%E2%80%93Fermi_model). Physics remembers him for his contribution, but he is most well-known for his algorithm, I'd argue (and I have not way of proving it, but its the 21st century, where facts and stats don't matter anymore it sems ...)

What was his algorithm, you ask? Well, he introduced the, drum roll please ... Gaussian elimination! Tada! Yes, Thomas developed a simplified version of the Gaussian elimination, where have I heard that before? It seems all direct methods are based on the Gaussian elimination, and people still get recognition for it.

To be fair, Thomas wrote down a simplification that applies to sparse matrices. Given that we do treat sparse matrices as their own sub-genre of matrices, I think it is fair to allow Thomas' algorithm as a standalone development. In fact, I insist we do, because this will ensure that we remember him for his simplification he made in a 1949 paper. CFD remembers him as one of the pioneers of CFD, when in reality, Thomas probably didn't have any idea what CFD was up until his death in 1992!

Thomas died in April 1992, to be precise, one month after the [Saab 2000](https://en.wikipedia.org/wiki/Saab_2000) had its maiden flight, and 7 months before the [Airbus 330](https://en.wikipedia.org/wiki/Airbus_A330) took to the skies. Seriously, in his life time, he witnessed the first powered flight and its development all the way to the A330. Why would you choose quantum physics? What a waste of time and talent ... 

So, let us look at his algorithm, then. There are two names in the literature you will find. Either people refer to it simply as the Thomas algorithm or the tri-diagonal matrix algorithm (TDMA). The T in TDMA gives us a hint of what type of matrices we are dealing with.

We looked at tri-diagonal matrices in this article already, but let's write out a simple system to see why it is called this way. Consider the following example:

$$
\begin{bmatrix}
D_1 & \beta_1 & 0 & \cdots & 0 \\
\alpha_2 & D_2 & \beta_2 & \cdots & 0 \\
0 & \alpha_3 & D_3 & \beta_3 & 0 \\
\vdots & \vdots & \vdots & \ddots & \beta_{n-1} \\
0 & 0 & 0 & \alpha_n & D_n
\end{bmatrix}
\begin{bmatrix}
\phi_1 \\
\phi_2 \\
\phi_3 \\
\vdots \\
\phi_n
\end{bmatrix}=
\begin{bmatrix}
b_1 \\
b_2 \\
b_3 \\
\vdots \\
b_n
\end{bmatrix}
\tag{eq:tdma-basic-form}
$$

This type of equation often arises in CFD. For example, if we consider the the 1D Poisson equation:

$$
\nabla^2 \phi = \mathbf{b}\\[1em]
\frac{\partial^2 \phi}{\partial x^2} = \mathbf{b} \\[1em]
\frac{\phi_{i+1}-2\phi_i + \phi_{i-1}}{(\Delta x)^2}=b_i \\[1em]
\left[\frac{1}{(\Delta x)^2}\right]\phi_{i+1} + \left[\frac{-2}{(\Delta x)^2}\right]\phi_{i} + \left[\frac{1}{(\Delta x)^2}\right]\phi_{i-1} = b_i
$$

In this discretisation, we have:

$$
\alpha = \left[\frac{1}{(\Delta x)^2}\right],\qquad D=\left[\frac{-2}{(\Delta x)^2}\right],\qquad \beta=\left[\frac{1}{(\Delta x)^2}\right]
$$

The tri-diagonal structure automatically appears for 1-dimensional problems, where our discretisation scheme uses 3 points to approximate derivatives (hence the name tri in TDMA). We have a diagonal and then two off-diagonals to either side. If this is given, then Thomas is arguing that we can simplify this form.

I already mentioned that we use the Gaussian elimination here, so how do we go about things? Well, you probably remember that we formed an upper triangular matrix, both in the Gaussian elimination and the LU decomposition. Since we have a tri-diagonal matrix, we can exploit this structure and get rid of elements in the lower triangular matrix.

Since we only have a single entry in the lower triangular matrix for each row, we only need to perform a single elimination step. For that, we define the mulitplier:

$$
m_i = \frac{\alpha_i}{D_{i-1}},\quad\text{for}\, i=2,...,n 
$$

We use this to modify our diagonal and right-hand side vector as:

$$
\tilde{D}_i = D_i - m_i \beta_{i-1}\\
\tilde{b}_i = b_i - m_i b_{i-1}
\tag{eq:tdma-substitutions}
$$

Let's write this out to see that this is indeed just Gaussian elimination and it produces an upper triangular matrix. Looking back at Eq.(\ref{eq:tdma-basic-form}), the first row can be written explicitly as:

$$
D_1 x_1 + \beta_1 x_2 = b_1
$$

The second row can be written as:

$$
\alpha_2 x_1 + D_2 x_2 + \beta x_3 = b_2
$$

Let's compute the multiplier $m_2$. This is:

$$
m_2=\frac{\alpha_i}{D_{i-1}}=\frac{\alpha_2}{D_{1}}
$$

We want to eliminate the lower triangular matrix, which in this case contains the value $\alpha_2$ in the second row. Thus, we take the second row and subtract from it the first row, multiplied by $m_2$. This gives us:

$$
(\alpha_2 x_1 + D_2 x_2 + \beta x_3) - \frac{\alpha_2}{D_{1}}\left(D_1 x_1 + \beta_1 x_2\right) = b_2 - \frac{\alpha_2}{D_{1}} b_1 \\[1em]
\alpha_2 x_1 - \frac{\alpha_2}{D_{1}}(D_1 x_1) + D_2 x_2 - \frac{\alpha_2}{D_{1}}\beta_1 x_2 + \beta x_3 = b_2 - \frac{\alpha_2}{D_{1}} b_1 \\[1em]
x_1\left(\alpha_2 - \frac{\alpha_2}{D_{1}}D_1\right) + x_2\left(D_2 - \frac{\alpha_2}{D_{1}}\beta_1\right) + \beta x_3 = b_2 - \frac{\alpha_2}{D_{1}} b_1 \\[1em]
x_1\left(\alpha_2 - \alpha_2\right) + x_2\left(D_2 - \frac{\alpha_2}{D_{1}}\beta_1\right) + \beta x_3 = b_2 - \frac{\alpha_2}{D_{1}} b_1 \\[1em]
x_1\left(0\right) + x_2\left(D_2 - \frac{\alpha_2}{D_{1}}\beta_1\right) + \beta x_3 = b_2 - \frac{\alpha_2}{D_{1}} b_1 \\[1em]
x_2\left(D_2 - \frac{\alpha_2}{D_{1}}\beta_1\right) + \beta x_3 = b_2 - \frac{\alpha_2}{D_{1}} b_1 \\[1em]
$$

Thus, we can see that by by multiplying the first equation by $m_i$ and subtracting it from the second, we do indeed get rid of the lower triangular matrix (the term containing $\alpha_2$, i.e. see Eq.(\ref{eq:tdma-basic-form})). We can see that the diagonal term is modified as:

$$
D_2 - \frac{\alpha_2}{D_{1}}\beta_1 \\[1em]
D_i - \frac{\alpha_i}{D_{i-1}}\beta_{i-1} \\[1em]
D_i - m_i\beta_{i-1} \\[1em]
$$

And, the right-hand side vector is modified as:

$$
b_2 - \frac{\alpha_2}{D_{1}} b_1 \\[1em]
b_i - \frac{\alpha_i}{D_{i-1}} b_{i-1} \\[1em]
b_i - m_i b_{i-1} \\[1em]
$$

This confirms that the relations given in Eq.(\ref{eq:tdma-substitutions}) are indeed correct and eliminate the lower triangular matrix. Since we only have a single off-diagonal to the left of the diagonal, we can use Eq.(\ref{eq:tdma-substitutions}) for all rows (starting at the second) row. This step is known as the forward elimination.

Once we have done that, we have transformed our original tri-diagonal matrix into an upper triangular matrix of the form:

$$
\begin{bmatrix}
\tilde{D}_1 & \beta_1 & 0 & \cdots & 0 \\
0 & \tilde{D}_2 & \beta_2 & \cdots & 0 \\
0 & 0 & \tilde{D}_3 & \beta_3 & 0 \\
\vdots & \vdots & \vdots & \ddots & \beta_{n-1} \\
0 & 0 & 0 & 0 & \tilde{D}_n
\end{bmatrix}
\begin{bmatrix}
\phi_1 \\
\phi_2 \\
\phi_3 \\
\vdots \\
\phi_n
\end{bmatrix}=
\begin{bmatrix}
\tilde{b}_1 \\
\tilde{b}_2 \\
\tilde{b}_3 \\
\vdots \\
\tilde{b}_n
\end{bmatrix}
$$

Doing a backward substituion again, we can write the last row of this modified matrix as:

$$
\tilde{D}_n \phi_n = \tilde{b}_n
$$

Once we know this value, we can compute any value in the rows above it as:

$$
\phi_i = \frac{\tilde{b}_i - \beta_i \phi_{i+1}}{\tilde{D}_i}
$$

In this way, we can loop from the back of the vector to the front, and obtain the values of $\phi_i$ in each loop.

Thomas' algorithm is extremely lucrative because it exploits the Gaussian elimination for sparse matrices, transforming the Gaussian elimination from an algorithm requiring $n^3$ operations down to $2n$. That is, we do one forward elimination step, where we get rid of all of the $\alpha_i$ values, which takes $n$ operations, and then we do one backward substitution step, which takes $n$ operations as well.

Keep in mind that this is a direct method; if we could, we would always use direct methods. They solve our linear system of equations, i.e. $\mathbf{Ax}=\mathbf{b}$, exactly. The Thomas algorithm, simply by exploiting the fact that our matrix $\mathbf{A}$ is tri-diagonal now, reduces the cost so much that there is no better algorithm available.

I could leave it here, finish the article, and pretend this is what we do in all of our CFD solvers, but life isn't that simple. What I have shown you here works for one dimensional problems. If you are dealing with one dimensional problems only, there likely isn't any other, better, faster, and more accurate algorithm out there!

The Thomas algorithm can be extended to two and three dimensions, but it becomes an iterative algorithm. To see this, let's write the Poisson example from before as a two dimensional problem. We have:

$$
\nabla^2 \phi = \mathbf{b}\\[1em]
\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = \mathbf{b} \\[1em]
\frac{\phi_{i+1,j}-2\phi_{i,j} + \phi_{i-1,j}}{(\Delta x)^2} + \frac{\phi_{i,j+1}-2\phi_{i,j} + \phi_{i,j-1}}{(\Delta y)^2}=b_i \\[1em]
\left[\frac{1}{(\Delta x)^2}\right]\phi_{i+1,j} + \left[\frac{-2}{(\Delta x)^2} + \frac{-2}{(\Delta y)^2}\right]\phi_{i,j} + \left[\frac{1}{(\Delta x)^2}\right]\phi_{i-1,j} + \left[\frac{1}{(\Delta y)^2}\right]\phi_{i,j+1} + \left[\frac{1}{(\Delta y)^2}\right]\phi_{i,j-1}= b_i
$$

When we build the coefficient matrix now, we have 5 entries per row, instead of three. So, instead of a tri-diagonal matrix, we have a penta-diagonal matrix, that is, a matrix with one main diagonal and 4 off diagonals. In practice, there is nothing wrong with extending the Thomas algorithm to penta-diagonal systems (i.e. 2D applications), however, this isn't done in practice.

What we would do instead is to cheat and write the above implicit system as a tri-diagonal matrix system isntead. This can be written as:

$$
\left[\frac{1}{(\Delta x)^2}\right]\phi_{i+1,j}^{k+1} + \left[\frac{-2}{(\Delta x)^2} + \frac{-2}{(\Delta y)^2}\right]\phi_{i,j}^{k+1} + \left[\frac{1}{(\Delta x)^2}\right]\phi_{i-1,j}^{k+1} = b_i - \left[\frac{1}{(\Delta y)^2}\right]\phi_{i,j+1}^k - \left[\frac{1}{(\Delta y)^2}\right]\phi_{i,j-1}^k
\tag{eq:penta-diagonal-implicit-x}
$$

or

$$
\left[\frac{1}{(\Delta y)^2}\right]\phi_{i,j+1}^{k+1} + \left[\frac{-2}{(\Delta x)^2} + \frac{-2}{(\Delta y)^2}\right]\phi_{i,j}^{k+1} + \left[\frac{1}{(\Delta y)^2}\right]\phi_{i,j-1}^{k+1}= b_i - \left[\frac{1}{(\Delta x)^2}\right]\phi_{i+1,j}^k - \left[\frac{1}{(\Delta x)^2}\right]\phi_{i-1,j}^k
\tag{eq:penta-diagonal-implicit-y}
$$

Here, I am indicating the values of $\phi^{k+1}$ are unknown, whereas values on the right-hand side at $\phi^k$ are known. We don't really know these values, but we can use them from the previous timestep or the initial conditions and we have a value for $\phi$ as an approximation. Then, we iterate over Eq.(\ref{eq:penta-diagonal-implicit-x}) or Eq.(\ref{eq:penta-diagonal-implicit-y}) until $\phi^{k+1}\approx\phi^k$ at each point $i,j$. 

Why is this beneficial? Well, we implement the Thomas algorithm once, and we can reuse it for two dimensional and three dimensional applications alike. We can also implement higher-order schemes, which have more than 3 points in their stencil and thus also destroy the tri-diagonal nature of out matrix.

While we are on this point, there is a method in CFD called the Alternating Direction Implicit (ADI) method. This method will first solve Eq.(\ref{eq:penta-diagonal-implicit-x}) and then Eq.(\ref{eq:penta-diagonal-implicit-x}), and repeat this until convergence, swapping between equations, which explains the alternating keyword, and the ADI method was (perhaps still is for some hardcore CFD developers) the goto choice to solve linear systems of equations.

So, we have established that the Gaussian elimination, while too expensive in its general form, can be reduced to more managable computational costs using either the LU decomposition or the Thomas algorithm. By combining Thomas with ADI, we have a solid framework for solving two and three dimensional equations. Is that the holy grail?

It depends. If you are happy to solve single-block, structured grid-type problems, then yes, to not read further. You will not need anything else, and just because the Thomas algorithm and ADI are dated methods, that doesn't make them anything less useful, or relevant. In fact, as already stated, for certain type of applications (mostly one dimensional problems), you can't beat Thomas in speed, storage requirements, accuracy, and general ignorance towards the aviation industry.

So, the way that I write about this method does suggest that, perhaps, we don't use it all that often in practice, and that is correct. And why don't we use Thomas anymore? Because we want to solve complex problems, which reuire an entire different data structure!

In general purpose CFD applications, we use unstructured grids, and these tend to destroy any diagonal characteristic of a matrix. While we are able to make a matrix look more like a diagonal matrix (using the [Cuthill McKee algorithm](https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm), for example), we cannot use the Thomas algorithm here as these matrices do not conform to the strict tri-diagonal structure we require for the Thomas algorithm.

Let me show you why, which is going to be much simpler than writing about it. The following sketch shows a very simple ustructured grid on the lect, where each cell as been given an ID. What I am showing on the right is the corresponding coefficient matrix that would arise as part of the discretisation, where an x denotes a non-zero entry. Where there are no entries, the values would be 0.

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\unstructured_grid_matrix_structure.png" alt="Sketch of a 2D domain on the left using unstructured triangular cells, with the corresponding coefficient matrix on the right showing non-zero entries" class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

We can see that each row corresponds to a cell, and we have a total of 8 cells. For each row (each cell), we look at which neighbour cells it has (which cells share a face with the cell we are currently focusing on). All neighbour cells will be part of the discretisation and contribute a non-zero term to the row. The column is equivalent to the cell ID of the neighbouring cells.

As we can see, the coefficient matrix that arises looks somewhat like a tri-diagonal system, but it clearly isn't. And, it only looks similar to a tri-diagonal system because the problem size if very small (8 cells). Make this problem larger, and you will see a lot of off-diagonal contributions.

So, We can't easily use Thomas and his algorithm here, and instead, we want to have algorithms that are robust and general, i.e. they can be used for one, two, and three dimensional flows, for single and multiblock, structured or unstructured applications alike.

We will start to develop these in the next section, but we will also see that their convergence properties aren't great. We will then try to formulate a different approach to solving these linear system of equation problems, which is where we will cross path again with our Krylov subspace, but I am getting ahead of myself. Let us first explore classical algorithm that can be applied to structured and unstructured grids alike!

## Matrix-free methods

So far, we have spent quite a bit of time preparing for solving lienar systems of equations, i.e. $\mathbf{Ax}=\mathbf{b}$. In this section, however, we will look how to solve this system without ever having to construct the coefficient matrix $\mathbf{A}$ in the first place. There are quite a few variants available, and we will look at the most commonly used and known algorithms here, as well as why we prefer some over others. 

### The Jacobi method 

The first method we will look at is the Jacobi method. It is so simple to implement that solving a linear system of equations really doesn't pose any difficulties. When you start to write your own solver and you have to solve an implicit system of equations, your usually start with the Jacobi method, because it is easy to implement and you are less likely to make any coding errors.

Having said that, the Jacobi method is also the slowest to convergence. However, we will see that the Jacobi method can be improved (the convergence rate can be doubled) by introducing a small fix. This small fix does not only double the convergence rate, but it also lowers storage requirements by a factor of 2. This is the Gauss-Seidel method we will look at in the next section.

I mention the Gauss-Seidel method here already, as it has a large overlap with the Jacobi method, just like the LU decomposition had a large overlap with the Gaussian elimination. And, while the Jacobi method is not used in practice (only, perhaps, as a starting point to implement other algorithms), the Gauss-Seidel method, which has a large overlap with the Jacobio method, is used in practice.

You will find it pretty much everywhere where a multigrid approach is used, which we will also discuss at the end of this section. Here, the Gauss-Seidel method is used as a so-called smoother. ANSYS Fluent computes its result with a Gauss-Seidel algorithm, OpenFOAM offers this option, too, and it is also a common option when we use multigrids in OpenFOAM.

You see, linear algebra without Gauss is quite impossible, and, presumably (a wild guess on my behalf), because Gauss' algorithms found so much widespread use in the early days of CFD, and its continous use in commercial CFD solvers, the Germans thought to honour Gauss by putting him on the 10 Deutsche Mark banknote, as seen below: 

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\gauss_money.jpg" alt="A 10 german mark banknote showing Gauss" class="wp-image-5550" style="width:600px"/><figcaption class="wp-element-caption">Figure reproduced from the<a href="https://www.bundesbank.de/en/tasks/cash-management/dem-banknotes-and-dem-coins/-/dem-banknotes-623628" target="_blank" rel="noopener" title="">Deutsche Bundesbank</a>.</figcaption></figure>
<!-- /wp:image -->

Just for reference, Germany had banknotes for 5, 10, 20, 50, 100, 200, 500, and 1000 (I think the 1000 banknote is a myth, I have never seen it). So, we put Gauss on the 10 banknote. Clearly, he was not as important as some of the others on more valuable banknotes (which were all poets, writers, and one architect, most of whom I have never heard of (there is a reason I don't live in Germany anymore ... I am not a model citizen ...)). Germany, *Land der Dichter und Denker*

Before we start to make our own currency now, and rank famous CFD people by importance, let's explore the core idea behind the Jacobi method. Let's say we have the following system of equations that we want to solve:

$$
\begin{bmatrix}
a_P & a_E & 0 & 0 & 0 \\[1em]
a_W & a_P & a_E & 0 & 0 \\[1em]
0 & a_W & a_P & a_E & 0 \\[1em]
0 & 0 & a_W & a_P & a_E \\[1em]
0 & 0 & 0 & a_W & a_P \\[1em]
\end{bmatrix}
\begin{bmatrix}
x_0 \\[1em]
x_1 \\[1em]
x_2 \\[1em]
x_3 \\[1em]
x_4 \\[1em]
\end{bmatrix}=
\begin{bmatrix}
b_0 \\[1em]
b_1 \\[1em]
b_2 \\[1em]
b_3 \\[1em]
b_4 \\[1em]
\end{bmatrix}
\tag{eq:jacobi-starting-problem}
$$

We saw this already in Eq.(\ref{eq:linear-system-general}), when we discretised a 1D equation. The Jacobi method also works for dense matrices, i.e. even though we use a sparse matrix here, we can use it just as well for any type of matrix (including those that arise from unstructured grid-based discretisations).

Let's look at the coefficient matrix $\mathbf{A}$ here first. What we do in the Jacobi method is to split this into its diagonal contribution, and its lower and upper triangular matrix. Let's do this first:

$$
\mathbf{A}=\mathbf{D}+\mathbf{L}+\mathbf{U}=\\[1em]
\underbrace{\begin{bmatrix}
a_P & a_E & 0 & 0 & 0 \\[1em]
a_W & a_P & a_E & 0 & 0 \\[1em]
0 & a_W & a_P & a_E & 0 \\[1em]
0 & 0 & a_W & a_P & a_E \\[1em]
0 & 0 & 0 & a_W & a_P \\[1em]
\end{bmatrix}}_\mathbf{A}=
\underbrace{\begin{bmatrix}
a_P & 0 & 0 & 0 & 0 \\[1em]
0 & a_P & 0 & 0 & 0 \\[1em]
0 & 0 & a_P & 0 & 0 \\[1em]
0 & 0 & 0 & a_P & 0 \\[1em]
0 & 0 & 0 & 0 & a_P \\[1em]
\end{bmatrix}}_\mathbf{D}+
\underbrace{\begin{bmatrix}
0 & 0 & 0 & 0 & 0 \\[1em]
a_W & 0 & 0 & 0 & 0 \\[1em]
0 & a_W & 0 & 0 & 0 \\[1em]
0 & 0 & a_W & 0 & 0 \\[1em]
0 & 0 & 0 & a_W & 0 \\[1em]
\end{bmatrix}}_\mathbf{L}+
\underbrace{\begin{bmatrix}
0 & a_E & 0 & 0 & 0 \\[1em]
0 & 0 & a_E & 0 & 0 \\[1em]
0 & 0 & 0 & a_E & 0 \\[1em]
0 & 0 & 0 & 0 & a_E \\[1em]
0 & 0 & 0 & 0 & 0 \\[1em]
\end{bmatrix}}_\mathbf{U}
\tag{eq:A-jacobi-decomposition}
$$

Now we use a trick: If you remember back to the beginning of this article, we looked at the main difference between explicit and implicit discretisations. We realised that explicit discretisations lead to a diagonal coefficient matrix, with only a single entry per row. Therefore, we were able to solve for the values in our unknown solution vector $\mathbf{x}$ (or $\mathbf{\phi}$) directly. When we looked at implicit systems of equations, we saw that we cannot do this, as we have now more than a single unknown per equation.

Using the matrix decomposition above, we can write this as:

$$
\mathbf{A}^{n+1}=\mathbf{D}^{n+1}+\mathbf{L}^{n+1}+\mathbf{U}^{n+1}
$$

This indicates that all of these coefficients came from terms that are unknwon, e.g. from $\phi_{i\pm 1,j\pm 1}^{n+1}$. In the Jacobi method, we now asssume that the lower and upper triangular matrix are known, which we can express as:

$$
\mathbf{A}^{n+1}=\mathbf{D}^{n+1}+\mathbf{L}^{n}+\mathbf{U}^{n}
\tag{eq:jacobi-core-idea}
$$

To make this clear, we can write this as:

$$
\mathbf{Ax}=\mathbf{b}\\[1em]
(\mathbf{D}+\mathbf{L}+\mathbf{U})\mathbf{x}=\mathbf{b}\\[1em]
\mathbf{D}\mathbf{x}^{n+1} + (\mathbf{L} + \mathbf{U})\mathbf{x}^{n}=\mathbf{b}
$$

If we write the equation in this way, we see that only the diagonal elements are unknown. Then, we can write the equation as:

$$
\mathbf{D}\mathbf{x}^{n+1}=\mathbf{b} - (\mathbf{L} + \mathbf{U})\mathbf{x}^{n}
\tag{eq:jacobi-core-idea-inserted}
$$

If we remember back at what we discussed in the beginning of this article, where we looked at the differences between explicit and implicit systems, we saw that explicit equations only have elements on the main diagonal, and so we can solve for them directly. Thus, using this decomposition of the Jacobi matrix into its diagonal and lower/upper triangular matrix, we have essentially just turned our implicit system into an explicit one. If we left things here, we would get a solution, but an explicit one.

This means we loose all favourable properties we gain with implicit systems, most notably that implicit system of equations tend to be unconditionally stable, i.e. regardless of the CFL number, we will get a solution. Explicit systems tend to be very sensitive to the CFL number and easily diverge if we don't set the timestep correctly based on the limitations imposed by an explicit discretisation. This can be shown by the [von Neumann stability analysis](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/explicit-vs-implicit-time-integration-and-the-cfl-condition/#aioseo-von-neumann), something which we have looked at previously.

So, we can't leave Eq.(\ref{eq:jacobi-core-idea-inserted}) as it is, but we can further refine this idea. Let us write this equation in a different form:

$$
\mathbf{D}\mathbf{x}^{k+1}=\mathbf{b} - (\mathbf{L} + \mathbf{U})\mathbf{x}^{k}
\tag{eq:jacobi-core-idea-solved-for-x}
$$

The only thing I have changed here is the superscript letter from $n$ to $k$. Since we want to solve for $\mathbf{x}^{k+1}$, let us also multiply by $\mathbf{D}^{-1}$ so that we isolate $\mathbf{x}^{k+1}$ on the left-hand side. This gives us:

$$
\mathbf{x}^{k+1}=\mathbf{D}^{-1}\left(\mathbf{b} - (\mathbf{L} + \mathbf{U})\mathbf{x}^{k}\right)
\tag{eq:jacobi-core-idea-solved-for-x-without-D}
$$

Instead of advancing the solution from $n$ to $n+1$, i.e. from the previous to the next time level, we take an iterative approach, where we start with $\mathbf{x}^{k}=\mathbf{x}^{n}$. Then, we can use this to compute $\mathbf{x}^{k+1}$. Now we can compare, is $\mathbf{x}^{k+1}\approx \mathbf{x}^{k}$? We have to specify how close the solution needs to be in order to be considered approximately equal (which we will do in a second).

If we find that $\mathbf{x}^{k+1}$ is quite far away from $\mathbf{x}^{k}$, well, then we set $\mathbf{x}^{k}=\mathbf{x}^{k+1}$, and we compute Eq.(\ref{eq:jacobi-core-idea-solved-for-x-without-D}) again. We do the same check again, i.e. is $\mathbf{x}^{k+1}\approx \mathbf{x}^{k}$? If it is, we stop and we set $\mathbf{x}^{n+1}=\mathbf{x}^{k+1}$, if it isn't, well, we set $\mathbf{x}^{k}=\mathbf{x}^{k+1}$ and evaluate Eq.(\ref{eq:jacobi-core-idea-solved-for-x-without-D}) again, until we have $\mathbf{x}^{k+1}\approx \mathbf{x}^{k}$.

This is an iterative procedure, which will force $\mathbf{x}^{k+1}$ and $\mathbf{x}^{k}$ to have the same value, at least to within some tolerance. If we go back to Eq.(\ref{eq:jacobi-core-idea-inserted}) now and we set $\mathbf{x}^{k+1}=\mathbf{x}^{n+1}$ and $\mathbf{x}^{k}=\mathbf{x}^{n}$, then we get:

$$
\mathbf{D}\mathbf{x}^{k+1}=\mathbf{b} - (\mathbf{L} + \mathbf{U})\mathbf{x}^{k}
$$

However, since $\mathbf{x}^{k+1}\approx \mathbf{x}^{k}$ after we have iterated for a sufficient number of iterations, we can write this also as:

$$
\mathbf{D}\mathbf{x}^{k+1}=\mathbf{b} - (\mathbf{L} + \mathbf{U})\mathbf{x}^{k+1}
$$

Or, using the equivalence $\mathbf{x}^{k+1}=\mathbf{x}^{n+1}$, we have:

$$
\mathbf{D}\mathbf{x}^{n+1}=\mathbf{b} - (\mathbf{L} + \mathbf{U})\mathbf{x}^{n+1}
$$

Thus, using an iterative approach, we have turned the otherwise explicit system (all variables on the right-hand side are evaluated at time level $n$) into a fully implicit system. This observation is key, and this will be used in the other methods in this subsection as well on matrix-free methods. We always have a system that looks like an explicit system, and we iterate until $\mathbf{x}^{k+1}\approx \mathbf{x}^{k}$.

To bring home this point, the following, somewhat pseudo C++ code shows how the Jacobi method could be implemented:

<!-- wp:kevinbatdorf/code-block-pro {"code":"for (int iteration = 0; iteration \u0026lt; maxJacobiIterations; ++iteration) {\n    // update solution\n    xOld = xNew;\n\n    // calculate new solution\n    for (int i = 0; i \u0026lt; numberOfPoints; ++i) {\n        // contributions from upper triangular matrix A \n        aE = A(i, i + 1);\n\n        // contributions from lower triangular matrix A \n        aW = A(i - 1, i);\n\n        // calculate new solution\n        xNew(i) = (b(i) - aE * xOld(i) - aW * xOld(i)) / A(i, i);\n    }\n    \n    // check convergence (e.g. eps = 10^-4)\n    if (L2Norm(xNew - xOld) \u0026lt; eps)\n        break;\n}","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e iteration = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; iteration \u0026lt; maxJacobiIterations; ++iteration) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // update solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    xOld = xNew;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; i \u0026lt; numberOfPoints; ++i) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from upper triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aE = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i + \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from lower triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aW = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i - \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003exNew\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) = (\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eb\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aE * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003exOld\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aW * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003exOld\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i)) / \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    }\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // check convergence (e.g. eps = 10^-4)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003eif\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eL2Norm\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(xNew - xOld) \u0026lt; eps)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003ebreak\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e}\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"cpp","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":20,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(2 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly>for (int iteration = 0; iteration &lt; maxJacobiIterations; ++iteration) {
    // update solution
    xOld = xNew;

    // calculate new solution
    for (int i = 0; i &lt; numberOfPoints; ++i) {
        // contributions from upper triangular matrix A 
        aE = A(i, i + 1);

        // contributions from lower triangular matrix A 
        aW = A(i - 1, i);

        // calculate new solution
        xNew(i) = (b(i) - aE * xOld(i) - aW * xOld(i)) / A(i, i);
    }
    
    // check convergence (e.g. eps = 10^-4)
    if (L2Norm(xNew - xOld) &lt; eps)
        break;
}</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> iteration = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; iteration &lt; maxJacobiIterations; ++iteration) {</span></span>
<span class="line"><span style="color: #6A9955">    // update solution</span></span>
<span class="line"><span style="color: #D4D4D4">    xOld = xNew;</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">    // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> i = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; i &lt; numberOfPoints; ++i) {</span></span>
<span class="line"><span style="color: #6A9955">        // contributions from upper triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aE = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i + </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // contributions from lower triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aW = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i - </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, i);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #DCDCAA">xNew</span><span style="color: #D4D4D4">(i) = (</span><span style="color: #DCDCAA">b</span><span style="color: #D4D4D4">(i) - aE * </span><span style="color: #DCDCAA">xOld</span><span style="color: #D4D4D4">(i) - aW * </span><span style="color: #DCDCAA">xOld</span><span style="color: #D4D4D4">(i)) / </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i);</span></span>
<span class="line"><span style="color: #D4D4D4">    }</span></span>
<span class="line"><span style="color: #D4D4D4">    </span></span>
<span class="line"><span style="color: #6A9955">    // check convergence (e.g. eps = 10^-4)</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">if</span><span style="color: #D4D4D4"> (</span><span style="color: #DCDCAA">L2Norm</span><span style="color: #D4D4D4">(xNew - xOld) &lt; eps)</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #C586C0">break</span><span style="color: #D4D4D4">;</span></span>
<span class="line"><span style="color: #D4D4D4">}</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

We see here that I have specified a variable called ```maxJacobiIterations```, so we need to tell our Jacobi method the maximum number of iterations we want to use. On line 3, we set the last computed solution ```xNew``` to the previous solution ```xOld```, which is the same as $\mathbf{x}^k=\mathbf{x}^{k+1}$ in the syntax we used before.

On line 6, we loop over all point in our mesh. If we have an unstructured grid, we typically only have a single loop, if we have a structured grid, we would have as many loops here as we have dimensions. In this case, I am using a 1D structured grid as an example, as we have done before.

On line 8 and 11, we get the coefficients from the upper and lower triangular matrix. Remember, these are comming directly from the coefficient matrix $\mathbf{A}$, as $\mathbf{A}=\mathbf{D}+\mathbf{L}+\mathbf{U}$ so instead of defining these lower and upper triangular matrices, I am just taking the values from $\mathbf{A}$ directly. See Eq.(\ref{eq:jacobi-starting-problem}) and Eq.(\ref{eq:A-jacobi-decomposition}) for the decomposition and where ```aE``` and ```aW``` are taken from.

On line 14, we compute a new solution for ```xNew```, or $\mathbf{x}^{k+1}$. This is the same as Eq.(\ref{eq:jacobi-core-idea-solved-for-x-without-D}). In this equation, we see that we have to multiply the right-hand side by $\mathbf{D}^{-1}$. Here, $\mathbf{D}$ is just the diagonal of the matrix. To invert the diagonal, we have $1/a_{ii}$, where $a_{ii}$ are the diagonal entries for each row $i$. We see that we do this by dividing the right-hand side by ```A(i, i)```, while the coefficients ```aE``` and ```aW``` would be coming from $\mathbf{L}$ and $\mathbf{U}$, respectively.

Once we have updated ```xNew```, i.e. $\mathbf{x}^{k+1}$, we go to line 18 and check the $||L_2||$ norm. If that is smaller than some small value ```eps```, typically in the range of $10^{-10}\le eps \le 10^{-3}$, then we say that ```xNew``` and ```xOld``` are approximately equal, or $\mathbf{x}^{k+1}\approx \mathbf{x}^{k}$. At this point, we break from the Jacobi iteration and continue with the rest of the simulation.

Just to be clear, we have said that the Jacobi method is computing the solution for

$$
\mathbf{x}^{k+1}=\mathbf{D}^{-1}\left(\mathbf{b} - (\mathbf{L} + \mathbf{U})\mathbf{x}^{k}\right)
$$

Evaluating this system until $\mathbf{x}^{k+1}\approx \mathbf{x}^{k}$ is equivalent to solving

$$
\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}
$$

Both equations will give us a solution to $\mathbf{x}$, but the Jacobi algorithm avoids having to invert the matrix $\mathbf{A}$ at the cost of having to iteratively find a solution to for $\mathbf{x}$. I just want to stress that even though these equations are different, they solve for the same quantity, so the Jacobi algorithm is just a reformulation of our original system $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$ into an easier system that we can solve.

So, we have seen that the Jacobi algorithm can be used to solve for our unknown vector $\mathbf{x}$ (or $\mathbf{\phi}$ in more general terms), and that we use iterations to get to the solution. But at what cost? There is an elegant way we can show mathematically how expensive the Jacobi iteration is. To do that, we have to introduce some new notation.

We saw that we can decompose our coefficient matrix as $\mathbf{A}=\mathbf{D}+\mathbf{L}+\mathbf{U}$. We then decomposed this further so that only the diagonal component is multiplied by $\mathbf{x}^{k+1}$, while both the lower and upper triangular matrix are multiplied by $\mathbf{x}^k$, and thus, these can be explicitly evaluated and put on the right-hand side of the equation. Let us write this out again. We have:

$$
\mathbf{Ax}=\mathbf{b}\\[1em]
(\mathbf{D}+\mathbf{L}+\mathbf{U})\mathbf{x}=\mathbf{b}\\[1em]
\mathbf{D}\mathbf{x}^{k+1} + (\mathbf{L}+\mathbf{U})\mathbf{x}^k=\mathbf{b}\\[1em]
\mathbf{D}\mathbf{x}^{k+1} =\mathbf{b} - (\mathbf{L}+\mathbf{U})\mathbf{x}^k
$$

In this form, we introduce now the following notation:

$$
\mathbf{D} = \mathbf{M},\qquad - (\mathbf{L}+\mathbf{U}) = \mathbf{N}
\tag{eq:M-N-matrices}
$$

Notice that the minus sign is part of the new matrix $\mathbf{N}$. The letter $\mathbf{M}$ and $\mathbf{N}$ may seem arbitrary, but these are used in the literature. The otherwise brilliant book by Moukalled *et al.* does not bother to specify what $\mathbf{M}$ and $\mathbf{N}$ is, though they do provide the decomposition of $\mathbf{A}=\mathbf{D}+\mathbf{L}+\mathbf{U}$. So, if you go through literature and find $\mathbf{M}$ and $\mathbf{N}$, then you know these are defined as in Eq.(\ref{eq:M-N-matrices}).

Let us insert this new definition, our linear system of equation becomes:

$$
\mathbf{M}\mathbf{x}^{k+1} =\mathbf{N}\mathbf{x}^k + \mathbf{b}
$$

Now, let's assume we have an exact equation (for example, from an analytic solution which may exist for a simple system). Let us write this as:

$$
\mathbf{M}\mathbf{x}^* =\mathbf{N}\mathbf{x}^* + \mathbf{b}
$$

Here, $\mathbf{x}^*$ states that the solution is exact (analytically known). In the previous two equations, the matrices $\mathbf{M}$ and $\mathbf{N}$ are the same. If we subtract the equation with the exact solution $\mathbf{x}^*$ from the equation with the approximate solution, then we get:

$$
\mathbf{M}\mathbf{x}^{k+1} - \mathbf{M}\mathbf{x}^* = \mathbf{N}\mathbf{x}^k - \mathbf{N}\mathbf{x}^* + \mathbf{b} - \mathbf{b}
\mathbf{M}(\mathbf{x}^{k+1} - \mathbf{x}^*) = \mathbf{N}(\mathbf{x}^k - \mathbf{x}^*)
\tag{eq:jacobi-with-M-and-N}
$$

Since $\mathbf{x}^*$ is the exact solution, the differences between $\mathbf{x}^{k+1}$ or $\mathbf{x}^{k}$ with $\mathbf{x}^*$ can be interpreted as an error, since we know the exact value. Let us, therefore, define:

$$
e^{k+1} = \mathbf{x}^{k+1} - \mathbf{x}^*\\[1em]
e^k = \mathbf{x}^k - \mathbf{x}^*
$$

We can insert this definition for the error into Eq.(\ref{eq:jacobi-with-M-and-N}) and get:

$$
\mathbf{M}e^{k+1} = \mathbf{N}e^k
$$

We can multiply the matrix $\mathbf{M}$ onto the right-hand side by forming its inverse and have:

$$
e^{k+1} = \mathbf{M}^{-1}\mathbf{N}e^k
$$

We further introduce the shorthand notation $\mathbf{G}=\mathbf{M}^{-1}\mathbf{N}$ and can simplify this equation as:

$$
e^{k+1} = \mathbf{G}e^k
\tag{eq:jacobi-error}
$$

So what does this equation tell us? Based on some properties of $\mathbf{G}$, we can make some judgement on how the error will gorw or shrink between iterations. The matrix $\mathbf{G}$ will have as many rows and columns as there are grid points in our simulation, and so we can compute the error increase/decrease for each row. Though, in reality, this isn't really practical. If we have a 1 million cell grid, then we probably don't want to evaluate 1 million errors individially.

Instead, we want some representative value for the matrix $\mathbf{G}$, and this is where the spectral radius, which have discussed before, comes in. As a short reminder, the spectral radius was given in Eq.(\ref{eq:spectral-radius}) as:

$$
\rho(\mathbf{A}) = \text{max}\{|\lambda_i|\}
$$

So, if we apply that to our new matrix $\mathbf{G}$, we have:

$$
\rho(\mathbf{G}) = \text{max}\{|\lambda_i|\}
$$

The spectral radius will look at all eigenvalues and return the largest (absolute) value. We can use this eigenvalue to state something about the convergence speed of the Jacobi method. For example, if the spectral radius of some matrix $\mathbf{G}$ is found to be $\rho(\mathbf{G}) = 0.5$, then, according to Eq.(\ref{eq:jacobi-error}), the error at the next iteration $k+1$ will be half of that at iteration $k$.

On the other hand, if $\rho(\mathbf{G}) = 0.1$, then the error at $k+1$ will reduce by a factor of 10 compared to the error at the previous iteration $k$. Thus, the smaller the spectral radius, the faster the error reduction. Remember, the error is a measure between the currently computed solution $\mathbf{x}^{k+1}$ and the analytically known (correct) solution $\mathbf{x}^*$. So, a reduction in error will bring $\mathbf{x}^{k+1}$ closer to $\mathbf{x}^*$.

As a result, once we get closer and closer to $\mathbf{x}^*$, the differences between $\mathbf{x}^{k+1}$ and $\mathbf{x}^{k}$ will become smaller and smaller, and so we reach a converged solution. This means that as long as $\rho(\mathbf{G})$ is smaller than 1, we get convergence, and the closer we are to 0, the faster we converge. However, we can also see that if $\rho(\mathbf{G})\gt 1$, the error will increase for each iteration, and we get divergence in our iterative solution procedure. We can write this in a compact form as:

$$
\rho(\mathbf{G})=
\begin{cases}
\lt 1 &\text{convergence}\\
1 &\text{stagnation}\\
\gt 1 &\text{divergence}
\end{cases}
$$

As it turns out, the Jacobi method is simple to implement, perhaps even easy to understand (that is something only you can judge), but this simplicity comes at a cost; it is just awefully slow! (as in, slow factorial ... OK, I stop with the nerdy math jokes ...). And, if we look at the definition for $\mathbf{G}$, we see that there isn't much we can do about it:

$$
\mathbf{G} = \mathbf{M}^{-1}\mathbf{N} = \mathbf{D}^{-1}(-(\mathbf{L}+\mathbf{U}))
$$

Since $\mathbf{G}$ depends entirely on the discretised system and the coefficient matrix that comes out of it, we can't modify this matrix and so the covnergence rate of the Jacobi method is fixed and can't be improved. Thus, if we wanted to improve convergence, we would have to change $\mathbf{G}$, and this is what this bloke on the 10 Deutsche Mark banknote did. Let's look at his idea next.

### The Gauss-Seidel method

Our story starts in the summer of 1821, in the Kingdom of Hannover. It was created along 38 other, sovereign, states within the German Confederation, in the wake of the collapse of the Holy Roman Empire (which existed for pretty much 1000 years, covering what now is known as Germany), all thanks to the Napoleonic Wars, which were fought by the French between 1803-1815. A small little French man was able to bring an Empire to its knees in a decade that existed for a millenial. I believe this is the reason *we all love* the French in Europe, right?

In 1819, a conference was held in Karlsbad by the German Confederation, and the Carlsbad Decrees were signed. This aimed at cracking down on nationalism and introducing stronger censorship to supress liberal ideas. Universities became under scrutiny and liberal professors were removed. Amids all of this, the government of the Kingdom of Hannover, which was geographically in what is now Germany but ruled by George IV (the King of the United Kingdom, *makes sense*), approached Gauss with a Task: They wanted to know how big their Kingdom was.

I'd love to think this was purely done for flexing, with George IV sending Twitter messages to his King friends asking *How big is yours?* (Kingdom, what did yo think?), but, with state lines shifting, the government needed to know what belonged to their Kingdom so they could tax people within the Kingdom correctly. Certainly, the military also had an interest to know where its borders were.

So, we had the Carlsbad Decrees, specifically targeting university professors (which Gauss was at that time), which, in many ways, was Germany's first attempt to establish a North Korean like regime (what is the difference between a dictator and a king with strong censorship?). We tried again some 100 years later with east Germany, proper communism and all of that Jazz, but we failed for a second time to establish a North Korea in the heart of Berlin.

As they say: [Third time's a charm](https://www.merriam-webster.com/dictionary/the%20third%20time%20is%20the%20charm), so who knows, perhaps there is a third, successful attempt in the year 2130 in whih Germany will have become a worthy successor of the North Korean *Empire* (at which point, we will be driving North Korean branded Hyundai and Kia cars). I really should stop now, I don't want to become the next [Jonny Somali](https://en.wikipedia.org/wiki/Johnny_Somali#South_Korea) ...

So, the scene is set, and Gauss began to measure the Kingdom. Similar efforts were completed earlier when the Kingdom of Denmark was measured and Gauss received help from them. At the same time, other states wanted to know their Kingdom size and boundaries to properly tax their citizens. But the way Gauss did it was special. Gauss didn't just measure the Kingdom of Hannover; he invented a whole new scientific disciplines.

By the time Gauss started his measurement, he was already famous like Einstein (of course, people would not compare him to Einstein back then, he wasn't born yet). The way he went about his measurements was to measure a triangle on a flat surface and use that as his baseline. Knowing the exact properties of one triangle laid out on a flat surface (length of each side and angles), Gauss could then determine properties of other triangles of much larger size.

So far, we have only used trigonometry, i.e. Gauss could stand on a tall building, or tower, and measure the angle between two points he could observe in the distance. He would repeat the same measurement on the other points, and then he could determine the area of each triangle. But, there is one problem. Let's look at his measurement in the following sketch:

<!-- wp:image {"width":"500px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\triangulation_measurement.png" alt="A sketch of a triangle where each angle is measured with some inaccuracies" class="wp-image-5550" style="width:500px"/></figure>
<!-- /wp:image -->

On the left, we see three points. Let's say we want to compute the area of the triangle that is spanned by these three points. If we follow Gauss, then we could look for three tall buildings, climb on top of each of them, and then measure the angle between each. Even today, we would make a small error in determining the angle between these three points, imagine how Gauss must have felt doing all by hand.

As we can see from the sketch, we will make a small error in determining the angle. In fact, we can get an idea for how accurate our measurement is by adding up all angle, i.e. $\alpha + \beta + \gamma$. Within a triangle, we know that all of that should add up to $180^\circ$, and it is unlikely that we get this value from our hand measurements.

Now, Gauss did not just span one triangle over the entire Kingdom of Hannover, but there were several triangles, leading to a network of triangles, as shown on the right. Here, all the angles between edges that meet at vertices will have a small error. We could leave it at that and accept that we will get a small error in our calculation, but Gauss wasn't happy with that.

Instead, he knew that all angle in a triangle should add up to $180^\circ$ and he used that to his advantage. But, since edges between triangles are shared, making changes to one triangle (making changes to one of the angles), will influence neighbouring triangles as well. So, the problem becomes one of minimising the error globably for all triangles when adjusting the angles.

Thankfully, while the Napoleonic Wars were in full swing, Gauss was busy publishing a method that we know as the least square problem nowadays. He could use the least square problem here as well to minimise measurement inaccuracies. We looked at the Least square problem already, when we looked at [Eigen and how to solve it with code](https://cfd.university/learn/essential-libraries-for-cfd-solver-development/eigen-the-swiss-arm-knife-of-cfd-libraries/#aioseo-computing-gradients-on-unstructured-grids-using-the-least-squares-and-related-approaches). You can have a quick look at it, but you won't need it to follow the rest of the discussion.

However, what I did show in the article, is that a least square problem can be formulated as:

$$
\mathbf{Ax}=\mathbf{b}
$$

Well, that is the starting point, but since the matrix $\mathbf{A}$ may not be a square matrix (i.e. the number of rows and column may not be the same), we have to make it square first. We do that by multiplying both sides by the transpose of $\mathbf{A}$ and have:

$$
\mathbf{A}^T\mathbf{Ax}=\mathbf{A}^T\mathbf{b}
$$

What we can see is that the solution to the least square problem can be formulated in a linear system of equation form $\mathbf{Ax}=\mathbf{b}$. And so, Gauss needed to a way to solve this equation. Even for him, Gaussian elimination was too expensive based on all of the triangles he had, and so he needed a cheaper solution. Krylov wasn't around, but Gauss, in a typical Gauss fashion, simply invented a new method to help him solve his problems, which he never published.

However, he used this method to correct his triangle network so to get a more accurate reading on the area of the Kingdom of Hannover. In fact, remember the banknote I showed you earlier? Well, if we turn it around, this is what we see:

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\gauss_money_back.jpg" alt="A 10 german mark banknote showing the backside of the banknote, including Gauss' triangulation of the Kingdom of Hannover" class="wp-image-5550" style="width:600px"/><figcaption class="wp-element-caption">Figure reproduced from the<a href="https://www.bundesbank.de/en/tasks/cash-management/dem-banknotes-and-dem-coins/-/dem-banknotes-623628" target="_blank" rel="noopener" title="">Deutsche Bundesbank</a>.</figcaption></figure>
<!-- /wp:image -->

The bottom right corner shows Gauss' triangular network of the Kingdom of Hannover. I only found out now. The Deutsche Mark was replaced by the Euro when I was 13 years old, so I suppose I did not have the strongest of interest in Gauss' work back then, so I never noticed. I did, however, read the book [Die Vermessung der Welt](https://www.amazon.co.uk/Die-Vermessung-Welt-Daniel-Kehlmann/dp/3499013223#), which I thought was quite good back then, still I had no idea about Gauss' triangulation (and the book literally deals with Gauss's attempt to measure the Kingdom).

So, what was that method I referred to? The method is known as the Gauss-Seidel method nowadays, and you will find it in general purpose CFD solvers. Gauss had no interest in publishing this method, he didn't see it as important. Some 50 years later, Seidel decided that perhaps some people may be interested in this method, and so he published it. In a true academic fashion, one person does all the work, and then the vultures come in, add their name to a paper, and get famous alongside. Academia isn't broken ...

In any case, let us have a look at the Gauss-Seidel method then. We can do that by looking at the starting point of the Jacobi method. We said that we can decompose a coefficient matrix $\mathbf{A}$ into its main diagonal, as well as its lower and upper triangle. We can use this decomposition to solve the linear system of equations $\mathbf{Ax}=\mathbf{b}$ as:

$$
\mathbf{Ax}=\mathbf{b}\\[1em]
(\mathbf{D}+\mathbf{L}+\mathbf{U})\mathbf{x}=\mathbf{b}\\[1em]
\mathbf{D}\mathbf{x}^{n+1} + (\mathbf{L} + \mathbf{U})\mathbf{x}^{n}=\mathbf{b}
\tag{eq:comparison-jacobi-gs-1}
$$

This was Jacobi. For the Gauss-Seidel algorithm, we make one small change:

$$
\mathbf{Ax}=\mathbf{b}\\[1em]
(\mathbf{D}+\mathbf{L}+\mathbf{U})\mathbf{x}=\mathbf{b}\\[1em]
(\mathbf{D} + \mathbf{L})\mathbf{x}^{n+1} + \mathbf{U}\mathbf{x}^{n}=\mathbf{b}
\tag{eq:comparison-jacobi-gs-2}
$$

Yes, that's it. This is the Gauss-Seidel algorithm in full display. I could stop here and move on, but it is worthwhile exploring why this method is better than the Jacobi method. First of all, we can look at the spectral radius. We said that we can decompose the coefficient matrix $\mathbf{A}$ into an $\mathbf{M}$ and $\mathbf{N}$ matrix, where $\mathbf{M}$ collects all matrices that multiply $\mathbf{x}^{n+1}$, while $\mathbf{N}$ collects all matrices that multiply $\mathbf{x}^n$. 

Before we do that, we have to put the known vector $\mathbf{x}^n$ on the right-hand side. Then, we can insert our definition for $\mathbf{M}$ and $\mathbf{N}$. This gives us:

$$
(\mathbf{D} + \mathbf{L})\mathbf{x}^{n+1} + \mathbf{U}\mathbf{x}^{n}=\mathbf{b}\\[1em]
(\mathbf{D} + \mathbf{L})\mathbf{x}^{n+1} = -\mathbf{U}\mathbf{x}^{n} + \mathbf{b}\\[1em]
\underbrace{(\mathbf{D} + \mathbf{L})}_{\mathbf{M}}\mathbf{x}^{n+1} = \underbrace{-\mathbf{U}}_\mathbf{N}\mathbf{x}^{n} + \mathbf{b}\\[1em]
\mathbf{M}\mathbf{x}^{n+1} = \mathbf{N}\mathbf{x}^{n} + \mathbf{b}
$$

Solving this for $\mathbf{x}^{n+1}$, we get:

$$
\mathbf{x}^{n+1} = \mathbf{M}^{-1}\mathbf{N}\mathbf{x}^{n} + \mathbf{M}^{-1}\mathbf{b}
$$

In the previous section on the Jacobi method, we saw that subtracting an exact equation from this equation produced an error equation of the form:

$$
e^{n+1} = \mathbf{M}^{-1}\mathbf{N}e^n\\[1em]
e^{n+1} = \mathbf{G}e^n
$$

If we look at the $\mathbf{G}$ equation, we have:

$$
\mathbf{G}_{GS} = \mathbf{M}^{-1}\mathbf{N} = (\mathbf{D} + \mathbf{L})^{-1}(-\mathbf{U})
$$

If we compare that with the $\mathbf{G}$ matrix we obtained fromt he Jacobi method, we had:

$$
\mathbf{G}_J = \mathbf{M}^{-1}\mathbf{N} = \mathbf{D}^{-1}(-(\mathbf{L}+\mathbf{U}))
$$

I am using the subscript $GS$ and $J$ to indicate that these are coming from the Gauss-Seidel and Jacobi method, respectively. So, while we can't really see if one of them is better than the other, we can say, however, that $\mathbf{G}_{GS} \ne \mathbf{G}_{J}$. And, therefore, we can also say that:

$$
\rho(\mathbf{G}_{GS}) \ne \rho(\mathbf{G}_{J})
$$

Well, but we still don't know which of these will give us a lower spectral radius. But we can build up a reasonable intuition by bring Eq.(\ref{eq:comparison-jacobi-gs-1}) and Eq.(\ref{eq:comparison-jacobi-gs-2}) into coefficient form. To get the coefficient form, we simply write out the explicit coefficients we use here. For example, the diagonal matrix $\mathbf{D}$ contains only elements on the main diagonal. If I wanted to index my matrix $\mathbf{A}$, and only get elements on the diagonal, well, I need to index rows and columns with the same index.

If my rows are indexed with the index $i$ and my columns are indexed with the index $j$, then, to get an element on the diagonal, $i$ must be equal to $j$, i.e. $i=j$. We could write this as $a_{ij}\forall i=j$, but that is very abstract, and we cannot really insert that into an equation. We can shortcut this and simply write $a_{ii}$. So, whatever row we are in, we are going to select the same column and are guaranteed to get the element on the diagonal.

So, we can start to transform Eq.(\ref{eq:comparison-jacobi-gs-1}) and Eq.(\ref{eq:comparison-jacobi-gs-2}) by removing the diagonal term and get:

$$
a_{ii}\mathbf{x}^{n+1} + (\mathbf{L} + \mathbf{U})\mathbf{x}^{n}=\mathbf{b}
(a_{ii} + \mathbf{L})\mathbf{x}^{n+1} + \mathbf{U}\mathbf{x}^{n}=\mathbf{b}
\tag{eq:from-matrix-to-coefficient-form-1}
$$

Now we need to think about the lower and upper matrix. To do that let's look back at Eq.(\ref{eq:A-jacobi-decomposition}), where we had:

$$
\mathbf{A}=\mathbf{D}+\mathbf{L}+\mathbf{U}=\\[1em]
\underbrace{\begin{bmatrix}
a_P & a_E & 0 & 0 & 0 \\[1em]
a_W & a_P & a_E & 0 & 0 \\[1em]
0 & a_W & a_P & a_E & 0 \\[1em]
0 & 0 & a_W & a_P & a_E \\[1em]
0 & 0 & 0 & a_W & a_P \\[1em]
\end{bmatrix}}_\mathbf{A}=
\underbrace{\begin{bmatrix}
a_P & 0 & 0 & 0 & 0 \\[1em]
0 & a_P & 0 & 0 & 0 \\[1em]
0 & 0 & a_P & 0 & 0 \\[1em]
0 & 0 & 0 & a_P & 0 \\[1em]
0 & 0 & 0 & 0 & a_P \\[1em]
\end{bmatrix}}_\mathbf{D}+
\underbrace{\begin{bmatrix}
0 & 0 & 0 & 0 & 0 \\[1em]
a_W & 0 & 0 & 0 & 0 \\[1em]
0 & a_W & 0 & 0 & 0 \\[1em]
0 & 0 & a_W & 0 & 0 \\[1em]
0 & 0 & 0 & a_W & 0 \\[1em]
\end{bmatrix}}_\mathbf{L}+
\underbrace{\begin{bmatrix}
0 & a_E & 0 & 0 & 0 \\[1em]
0 & 0 & a_E & 0 & 0 \\[1em]
0 & 0 & 0 & a_E & 0 \\[1em]
0 & 0 & 0 & 0 & a_E \\[1em]
0 & 0 & 0 & 0 & 0 \\[1em]
\end{bmatrix}}_\mathbf{U}
$$

Let's say I am in the third row ($i=3$). If I wanted to get the diagonal entry, we said we have $a_{ii}=a_{33}=a_P$. And if I wanted to get all elements in the lower triangular matrix? Well, we see that we have one entry to the left of $a_{33}$. That is, in the lower triangular matrix, we have $a_{32}=a_W$. In index form, we could write this as $a_{i,i-1}$.

Since the lower triangular matrix is, well, a matrix, we have to evaluate the matrix vector product. In Eq.(\ref{eq:comparison-jacobi-gs-2}), for example, we have to evaluate $\mathbf{L}\mathbf{x}^{n+1}$. If we are in row three in our example matrix above, then we can write this out explicitly as:

$$
\mathbf{L}\mathbf{x}^{n+1} = l_{31}x_1^{n+1} + l_{32}x_2^{n+1} = 0\cdot x_1^{n+1} + a_W\cdot x_2^{n+1} = a_W\cdot x_2^{n+1}
$$

This is now one example, but we want to generalise this. First, we notice that the column index of the matrix is the same as the index of our unknown vector $\mathbf{x}$. Next, we are in row three, so $i=3$. We also notice that we have a summation of two terms, so we want to express that with a summation. Let's give this a try. The above equation can be generalised as:

$$
\mathbf{L}\mathbf{x}^{n+1} = \sum_{j=1}^{i-1}l_{ij}x_j^{n+1}
$$

We start the summation in the first colum, hence the summation operation starts at $j=1$. We are in the third row, and the lower triangular matrix does not include the main diagonal, so the summation has to stop before the diagonal. The diagonal is at $i=3$, so the summation needs to stop at $i-1=2$ (and we saw that we only had two terms in our summation before, this part checks out).

Then, $l_{ij}$ are simply the coefficients in the lower triangular matrix, while $x_j$ are the components in our unknown solution vector $\mathbf{x}$. So, inserting this into Eq.(\ref{eq:from-matrix-to-coefficient-form-1}), we have:

$$
a_{ii}\mathbf{x}^{n+1} + \sum_{j=1}^{i-1}l_{ij}x_j^{n} + \mathbf{U}\mathbf{x}^{n}=\mathbf{b}
a_{ii}\mathbf{x}^{n+1} + \sum_{j=1}^{i-1}l_{ij}x_j^{n+1} + \mathbf{U}\mathbf{x}^{n}=\mathbf{b}
\tag{eq:from-matrix-to-coefficient-form-2}
$$

What about the matrix $\mathbf{U}$? Well, this is just the upper triangular matrix, so we have to do something very similar to the lower triangular matrix. Instead of summing from $j=1$ to $i-1$, which covers all entries in the lower triangular matrix up until the diagonal, we have to do the sum from $i+1$ to $N$, where $N$ is the number of rows and columns in our coefficient matrix $\mathbf{A}$. If we do that, we can write, for example:

$$
\mathbf{U}\mathbf{x}^{n} = \sum_{j=i+1}^{N}u_{ij}x_j^{n}
$$

The summation we have to evaluate hasn't changed, only the start and end of the summation. So, we can insert that into Eq.(\ref{eq:from-matrix-to-coefficient-form-2}) and obtain:

$$
a_{ii}\mathbf{x}^{n+1} + \sum_{j=1}^{i-1}l_{ij}x_j^{n}   + \sum_{j=i+1}^{N}u_{ij}x_j^{n}=\mathbf{b}
a_{ii}\mathbf{x}^{n+1} + \sum_{j=1}^{i-1}l_{ij}x_j^{n+1} + \sum_{j=i+1}^{N}u_{ij}x_j^{n}=\mathbf{b}
\tag{eq:from-matrix-to-coefficient-form-3}
$$

The only vector left in this expression is the right-hand side vector $\mathbf{b}$ and the unknown vector $\mathbf{x}^{n+1}$ multiplying the diagonal coefficient. Since we go through this equation row by row, i.e. we go through each row in our matrix with index $i$, we have to replace both $\mathbf{b}$ and $\mathbf{x}^{n+1}$ by their $i$-th component. Thus, we replace these vectors with $b_i$ and $x_i^{n+1}$, respectively, and we get:

$$
a_{ii}x_i^{n+1} + \sum_{j=1}^{i-1}l_{ij}x_j^{n}   + \sum_{j=i+1}^{N}u_{ij}x_j^{n}=b_i
a_{ii}x_i^{n+1} + \sum_{j=1}^{i-1}l_{ij}x_j^{n+1} + \sum_{j=i+1}^{N}u_{ij}x_j^{n}=b_i
\tag{eq:from-matrix-to-coefficient-form-4}
$$

Remember, the first equation is our Jacobi method, the second equation is our Gauss-Seidel equation. We can solve both of these now for $x_i^{n+1}$ and see what the resulting update equation is. I will split both equations now into their own equation so we can look at them separately. For the Jacobi method, we get:

$$
a_{ii}x_i^{n+1} + \sum_{j=1}^{i-1}l_{ij}x_j^{n}   + \sum_{j=i+1}^{N}u_{ij}x_j^{n}=b_i\\[1em]
x_i^{n+1} =\frac{1}{a_{ii}}\left(b_i - \sum_{j=1}^{i-1}l_{ij}x_j^{n} - \sum_{j=i+1}^{N}u_{ij}x_j^{n}\right)
\tag{eq:jacobi-coefficient-form}
$$

For the Gauss-Seidel method, we get:

$$
x_i^{n+1} = \frac{1}{a_{ii}}\left(b_i - \sum_{j=1}^{i-1}l_{ij}x_j^{n+1} - \sum_{j=i+1}^{N}u_{ij}x_j^{n}\right)
\tag{eq:gauss-seidel-coefficient-form}
$$

Hang on, for the Gauss-Seidel method, we now have $x_i^{n+1}$ on the left hand side and $x_j^{n+1}$ on the right-hand side. Didn't we just say that we wanted to collect all unknowns on the left-hand side? Yes, we did that. But why is $x_j^{n+1}$ on the right-hand side? Well, even though it is, at some point, unknown, once we evaluate $x_i^{n+1}$, it is already known.

Remember, $x_j^{n+1}$ multiplies only the lower triangular matrix, which goes from $j=1$ to $i-1$. This, when we evaluate $x_i^{n+1}$, say, for example, we are in the third row and have $i=3$, then we want to compute $x_3^{n+1}$. But, the lower triangular matrix goes only up until $i-1=3-1=2$. So, the lower matrix coefficients are multiplied by $x_1$ and $x_2$.

If we iterate over all rows in the matrix, typically starting at the beginning, i.e. row 1, and looping to the last row, then, at the time we evaluate $x_3$, we have already computed $x_1$ and $x_2$ in a previous iteration. Thus, once we compute $x_3^{n+1}$, we have already computed $x_1^{n+1}$ and $x_2^{n+1}$. Since we only need these previously computed values in our lower triangular coefficient matrix, we can write this on our right-hand side and assume these values are know.

This is the reason why the Gauss-Seidel method is faster, by approximately a factor of two, compared to the Jacobi method, because it uses the updated solution immediately in each iteration. The best part about the Gauss-Seidel method is that once you have implemented the Jacobi method, getting the Gauss-Seidel method implmeneted is straightforward.

Let's modify the pseudo C++ code we have seen before for the Jacobi method. To modify it and make it the Gauss-Seidel method, this is what we have to do: 

<!-- wp:kevinbatdorf/code-block-pro {"code":"for (int iteration = 0; iteration \u0026lt; maxJacobiIterations; ++iteration) {\n    // update solution\n    xOld = xNew;\n\n    // calculate new solution\n    for (int i = 0; i \u0026lt; numberOfPoints; ++i) {\n        // contributions from upper triangular matrix A \n        aE = A(i, i + 1);\n\n        // contributions from lower triangular matrix A \n        aW = A(i - 1, i);\n\n        // calculate new solution\n        xNew(i) = (b(i) - aE * xNew(i) - aW * xNew(i)) / A(i, i);\n    }\n    \n    // check convergence (e.g. eps = 10^-4)\n    if (L2Norm(xNew - xOld) \u0026lt; eps)\n        break;\n}","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e iteration = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; iteration \u0026lt; maxJacobiIterations; ++iteration) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // update solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    xOld = xNew;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; i \u0026lt; numberOfPoints; ++i) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from upper triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aE = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i + \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from lower triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aW = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i - \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003exNew\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) = (\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eb\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aE * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003exNew\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aW * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003exNew\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i)) / \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    }\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // check convergence (e.g. eps = 10^-4)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003eif\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eL2Norm\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(xNew - xOld) \u0026lt; eps)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003ebreak\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e}\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"cpp","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":20,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(2 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly>for (int iteration = 0; iteration &lt; maxJacobiIterations; ++iteration) {
    // update solution
    xOld = xNew;

    // calculate new solution
    for (int i = 0; i &lt; numberOfPoints; ++i) {
        // contributions from upper triangular matrix A 
        aE = A(i, i + 1);

        // contributions from lower triangular matrix A 
        aW = A(i - 1, i);

        // calculate new solution
        xNew(i) = (b(i) - aE * xNew(i) - aW * xNew(i)) / A(i, i);
    }
    
    // check convergence (e.g. eps = 10^-4)
    if (L2Norm(xNew - xOld) &lt; eps)
        break;
}</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> iteration = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; iteration &lt; maxJacobiIterations; ++iteration) {</span></span>
<span class="line"><span style="color: #6A9955">    // update solution</span></span>
<span class="line"><span style="color: #D4D4D4">    xOld = xNew;</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">    // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> i = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; i &lt; numberOfPoints; ++i) {</span></span>
<span class="line"><span style="color: #6A9955">        // contributions from upper triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aE = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i + </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // contributions from lower triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aW = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i - </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, i);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #DCDCAA">xNew</span><span style="color: #D4D4D4">(i) = (</span><span style="color: #DCDCAA">b</span><span style="color: #D4D4D4">(i) - aE * </span><span style="color: #DCDCAA">xNew</span><span style="color: #D4D4D4">(i) - aW * </span><span style="color: #DCDCAA">xNew</span><span style="color: #D4D4D4">(i)) / </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i);</span></span>
<span class="line"><span style="color: #D4D4D4">    }</span></span>
<span class="line"><span style="color: #D4D4D4">    </span></span>
<span class="line"><span style="color: #6A9955">    // check convergence (e.g. eps = 10^-4)</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">if</span><span style="color: #D4D4D4"> (</span><span style="color: #DCDCAA">L2Norm</span><span style="color: #D4D4D4">(xNew - xOld) &lt; eps)</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #C586C0">break</span><span style="color: #D4D4D4">;</span></span>
<span class="line"><span style="color: #D4D4D4">}</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

Did you spot the difference? The only change we have to do is on line 14, where we have replaced all ```xOld``` arrays by ```xNew``` arrays. Essentially, we are reading and writing here to the same array and use an updated solution here if already available. If not, well, then we simply take the solution that we had in the array from the previous iteration (or initial condition).

I do use both ```xNew``` and ```xOld``` here, but I only keep ```xOld``` around so that I can compute the difference between two iterations to check if the solution has converged. However, if we don't need to check for convergence, and we know we are only ever going to use it for a few iterations (for example, during a multigrid cycle), then we can remove the convergence checking and simplify the Gauss-Seidel method to:

<!-- wp:kevinbatdorf/code-block-pro {"code":"for (int iteration = 0; iteration \u0026lt; maxJacobiIterations; ++iteration) {\n    // calculate new solution\n    for (int i = 0; i \u0026lt; numberOfPoints; ++i) {\n        // contributions from upper triangular matrix A \n        aE = A(i, i + 1);\n\n        // contributions from lower triangular matrix A \n        aW = A(i - 1, i);\n\n        // calculate new solution\n        x(i) = (b(i) - aE * x(i) - aW * x(i)) / A(i, i);\n    }\n}","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e iteration = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; iteration \u0026lt; maxJacobiIterations; ++iteration) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; i \u0026lt; numberOfPoints; ++i) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from upper triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aE = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i + \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from lower triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aW = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i - \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) = (\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eb\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aE * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aW * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i)) / \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    }\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e}\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"cpp","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":13,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(2 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly>for (int iteration = 0; iteration &lt; maxJacobiIterations; ++iteration) {
    // calculate new solution
    for (int i = 0; i &lt; numberOfPoints; ++i) {
        // contributions from upper triangular matrix A 
        aE = A(i, i + 1);

        // contributions from lower triangular matrix A 
        aW = A(i - 1, i);

        // calculate new solution
        x(i) = (b(i) - aE * x(i) - aW * x(i)) / A(i, i);
    }
}</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> iteration = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; iteration &lt; maxJacobiIterations; ++iteration) {</span></span>
<span class="line"><span style="color: #6A9955">    // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> i = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; i &lt; numberOfPoints; ++i) {</span></span>
<span class="line"><span style="color: #6A9955">        // contributions from upper triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aE = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i + </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // contributions from lower triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aW = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i - </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, i);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i) = (</span><span style="color: #DCDCAA">b</span><span style="color: #D4D4D4">(i) - aE * </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i) - aW * </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i)) / </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i);</span></span>
<span class="line"><span style="color: #D4D4D4">    }</span></span>
<span class="line"><span style="color: #D4D4D4">}</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

In this case, we have only a single ```x``` array, which we constantly write to and read from at the same time. We have reduced the stroage requirements now by a factor of two (we only store one array, not two arrays), and we have gained convergence acceleration by a factor of about two. Not bad fort such a small change!

In the Gauss-Seidel algorithm, for example, in the pseudo code given above, we always looped from the first row of the coefficient matrix to the last row, or, in terms of our loop, we started at $i=0$ and went the ```numberOfPoints```. This is known as the **forward Gauss Seidel**. Now, here is a tough nut for you to crack. Can you imagine what the **backward Gauss Seidel** method is doing?

If you guessed going from the back (last row) to the front (first row), wow, you are on your way to a noble prize! (I mean, Gauss set ```xOld=xNew``` and called this a new method, Seidel didn't do anything beside publish this and we remember his name for it. If this is the bar for claiming academic fame, a noble prize can't be that difficult to obtain I recon ...)

Ok, we know what the forward and backward Gauss-Seidel method (forward and backward looping), but riddle me this: What is the **symmetric Gauss-Seidel** method? Well, this is a good application of [Occam's razor](https://en.wikipedia.org/wiki/Occam's_razor), which essentially states that the most logical solution is the simplest one. This is one example where we are allowed to use business mathematics:

$$
\text{symmetric Gauss-Seidel} = \text{forward Gauss-Seidel} + \text{backward Gauss-Seidel}
$$

So, the symmetric Gauss-Seidel method does a forward sweep first, followed by a backward sweep. Why is this a good idea? If we are only doing a forward or backward sweep, information will always propagate in one direction. However, if we do a symmetric sweep, information can propagate in different directions. To make this clear, take a look at the discretisation of the laplacian operator in one dimension, for example:

$$
\frac{\partial^2\phi}{\partial x^2}\approx \frac{\phi_{i+1}-2\phi_i + \phi_{i-1}}{(\Delta x)^2}=\phi_{i+1}\underbrace{\left[\frac{1}{(\delta x)^2}\right]}_\mathbf{U} + \phi_i\underbrace{\left[\frac{-2}{(\Delta x)^2}\right]}_\mathbf{D} + \phi_{i-1}\underbrace{\left[\frac{1}{(\Delta x)^2}\right]}_\mathbf{L}
\tag{eq:laplacian-example-sgs}
$$

Here, we can see that if we used a forward Gauss-Seidel method, for example, we would always update the solution such that $\phi_{i-1}$ would already be known (always, as this is coming from the lower triangular matrix) when we update $\phi_i$, so the solution for $\phi_{i+1}$ is always lagging behind. This asymmetry can then introduce a directional dependence which may slow down convergence.

The cure is to use a symmetric Gauss-Seidel method in this case, which will update $\phi_{i+1}$ and $\phi_{i-1}$ in an alternating fashion, aiding to preserve the character of the Laplacian operator (diffusion in each direction in equal parts, rather than stronger diffusion in one direction). These type of operators lead to symmetric matrices. 

We can see from Eq.(\ref{eq:laplacian-example-sgs}) that the contributions to the lower and upper triangular matrices $\mathbf{L}$ and $\mathbf{U}$ are going to be the same, so the coefficient matrix $\mathbf{A}$, which consists of the diagonal, upper, and lower triangular matrix, i.e. $\mathbf{A} = \mathbf{D} + \mathbf{L} + \mathbf{U}$, is going to be a symmetric matrix. For these type of matrices, the symmetric Gauss-Seidel method has computational advantages (it may converge faster) over the pure forward or backward Gauss Seidel method.

In pseudo C++ code we can write this as:

<!-- wp:kevinbatdorf/code-block-pro {"code":"for (int iteration = 0; iteration \u0026lt; maxJacobiIterations; ++iteration) {\n    // forward sweep\n    for (int i = 0; i \u0026lt; numberOfPoints; ++i) {\n        // contributions from upper triangular matrix A \n        aE = A(i, i + 1);\n\n        // contributions from lower triangular matrix A \n        aW = A(i - 1, i);\n\n        // calculate new solution\n        x(i) = (b(i) - aE * x(i) - aW * x(i)) / A(i, i);\n    }\n    \n    // backward sweep\n    for (int i = numberOfPoints - 1; i \u003e= 0; \u002d\u002di) {\n        // contributions from upper triangular matrix A \n        aE = A(i, i + 1);\n\n        // contributions from lower triangular matrix A \n        aW = A(i - 1, i);\n\n        // calculate new solution\n        x(i) = (b(i) - aE * x(i) - aW * x(i)) / A(i, i);\n    }\n}","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e iteration = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; iteration \u0026lt; maxJacobiIterations; ++iteration) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // forward sweep\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; i \u0026lt; numberOfPoints; ++i) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from upper triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aE = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i + \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from lower triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aW = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i - \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) = (\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eb\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aE * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aW * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i)) / \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    }\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e    // backward sweep\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e (\u003c/span\u003e\u003cspan style=\u0022color: #569CD6\u0022\u003eint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i = numberOfPoints - \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; i \u0026gt;= \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e; \u002d\u002di) {\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from upper triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aE = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i + \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // contributions from lower triangular matrix A \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        aW = \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i - \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e        // calculate new solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) = (\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eb\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aE * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i) - aW * \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003ex\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i)) / \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eA\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(i, i);\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    }\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e}\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"cpp","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":25,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(2 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly>for (int iteration = 0; iteration &lt; maxJacobiIterations; ++iteration) {
    // forward sweep
    for (int i = 0; i &lt; numberOfPoints; ++i) {
        // contributions from upper triangular matrix A 
        aE = A(i, i + 1);

        // contributions from lower triangular matrix A 
        aW = A(i - 1, i);

        // calculate new solution
        x(i) = (b(i) - aE * x(i) - aW * x(i)) / A(i, i);
    }
    
    // backward sweep
    for (int i = numberOfPoints - 1; i >= 0; --i) {
        // contributions from upper triangular matrix A 
        aE = A(i, i + 1);

        // contributions from lower triangular matrix A 
        aW = A(i - 1, i);

        // calculate new solution
        x(i) = (b(i) - aE * x(i) - aW * x(i)) / A(i, i);
    }
}</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> iteration = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; iteration &lt; maxJacobiIterations; ++iteration) {</span></span>
<span class="line"><span style="color: #6A9955">    // forward sweep</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> i = </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; i &lt; numberOfPoints; ++i) {</span></span>
<span class="line"><span style="color: #6A9955">        // contributions from upper triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aE = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i + </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // contributions from lower triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aW = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i - </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, i);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i) = (</span><span style="color: #DCDCAA">b</span><span style="color: #D4D4D4">(i) - aE * </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i) - aW * </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i)) / </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i);</span></span>
<span class="line"><span style="color: #D4D4D4">    }</span></span>
<span class="line"><span style="color: #D4D4D4">    </span></span>
<span class="line"><span style="color: #6A9955">    // backward sweep</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> (</span><span style="color: #569CD6">int</span><span style="color: #D4D4D4"> i = numberOfPoints - </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">; i &gt;= </span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">; --i) {</span></span>
<span class="line"><span style="color: #6A9955">        // contributions from upper triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aE = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i + </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // contributions from lower triangular matrix A </span></span>
<span class="line"><span style="color: #D4D4D4">        aW = </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i - </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, i);</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955">        // calculate new solution</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i) = (</span><span style="color: #DCDCAA">b</span><span style="color: #D4D4D4">(i) - aE * </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i) - aW * </span><span style="color: #DCDCAA">x</span><span style="color: #D4D4D4">(i)) / </span><span style="color: #DCDCAA">A</span><span style="color: #D4D4D4">(i, i);</span></span>
<span class="line"><span style="color: #D4D4D4">    }</span></span>
<span class="line"><span style="color: #D4D4D4">}</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

Finally, you may ask, what is business mathematics? Its when we add words together and pretend its maths. Next time you see someone use business math, you have my permission to write an integral infront of each word and integrate the constant 1 from 0 to 1 (which is, well, 1). We haven't changed anything, the equation is still correct, but business math people won't realise. So much much for business *math* ...

### Successive Over-relaxation (SOR)

Let's talk about under-relaxation for a moment. When we solve linear partial differential equations, we can simply form a discretised version of them and then compute an updated solution. Job done! But, with non-linear partial differential equations, small inaccuracies can quickly multiply and lead to divergence in our simulations.

If we look at the [SIMPLE](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-solve-incompressible-and-compressible-flows-in-cfd/#aioseo-pressure-projection-based-simple) method, which was introduced to solve incompressible flows and, in my view, is a blatant rip-off and plagiarised (and somehoe inferior) method based on Chorin's exact pressure projection method, Patankar and Spalding butchered Chorin's pressure projection method so badly that they would no longer give us an accurate prediction of the flow.

Some may say this was required to take away some of that non-linearity so as to stabilise the non-linear partial differential equations and, as a result, the flow field as well. However, I'd argue that we could have done the same with Chorin's exact pressure projection method but well, I don't teach at Imperial College so *obviously* my opinions aren't as valuable.

I mean, if cook and cut yourself in the finger, if you wanted to stop the bleeding, I'd argue applying a bandaid is a sensible solution. However, if we apply the same logic as ther SIMPLE algorithm here, we would have to chop off our hand. Sure, it will stop the bleeding in our finger, but we have just created a dumpster fire of other issues elsewhere.

So, the SIMPLE algorithm by itself doesn't work. We linearise an otherwise non-linear partial differential equation, and we even get rid of parts of the non-linear term, just to make our life easier (or *simpler*?). But then, we realise, we probably can do better, and we bring back these non-linear terms we have dropped, and this allows us to introduce a new method, called the SIMPLEC method, where the C stands for consistent. In reality, the SIMPLEC method only brings the SIMPLE method closer to Chorin's exact projection method.

In any case, let's focus in on the SIMPLE algorithm. It has so many simplification that we are only getting an approximate prediction of the flow (once we get convergence, we will recover a sensible flow field again, as it still retains Chorin's pressure projection core ideas). This means that if we were to update our solution with the SIMPLE algorithm and then continue to update our solution with these new updates fields, we would eventually get divergence in our simulations.

So, we introduce a fix, called under-relaxation. What we do here is to only update our new solution $\mathbf{x}^{n+1}$ by a fraction (somewhere between 0 and 1, or 0% and 100%). We write this as:

$$
\mathbf{x}^{n+1}=\alpha\mathbf{x}^{n+1} + (1-\alpha)\mathbf{x}^{1}
\tag{eq:under-relaxation}
$$

Here, $\alpha$ is known as an under-relaxation parameter, and it is taken values between $0\lt\alpha\le 1$. While Chorin's pressure projection method can operate with $\alpha=1$ pretty confidentally (eventually we may need to use some form of under-relaxation if we start to use turbulence models, which themselves are not exact), the SIMPLE algorithm requires heavey under-relaxation. Values of $\alpha_U=0.7$ and $\alpha_p=0.3$ are common for the velocity field and pressure field, respectively.

If we insert $\alpha=0.7$, for example, into Eq.(\ref{eq:under-relaxation}), we see that we are now taking only 0.7 of the newly computed solution, and we add $1-\alpha=1-0.7=0.3$ from the previously computed solution, i.e. the solution at time level $n$. So, we can see under-relaxation also as a weighted sum, or linear blend between two quantities.

This under-relaxation stabilises the SIMPLE algorithm suifficiently so that solutions can be obtained again. As the name *under-relaxation* suggest, values of $\alpha$ are in the range of $0\lt\alpha\le 1$, as already highlighted. But, we can also introduce over-relaxation, where we require $\alpha\gt 1$. If we tried to do that, we would get divergence, and so, over-relaxation isn't something we commonly use in CFD. Except for one case: The Gauss-Seidel method.

Yes, this is now, what, the 5th section or so where we talk about Gauss again? I told you, we can't do linear algebra without Gauss, just like we can't have the UK without mildly offensive boarder security guards at heathrow airport. They *love* welcoming foreigners (with fully settled status) into their little Kingdom ... It is just part of the deal (and I shall avoid heathrow airport for a while). Is now the right time to talk about the time I smuggled German sausages into Korea? Perhaps a story for another time ...

So, let us apply the same idea of under-relaxation to the Gauss-Seidel method. I will then show you how this under-relaxation can actually lead to stable results with over-relaxation. Happy days!

Let's look at the Gauss-Seidel algorithm again, which was given in Eq.(\ref{eq:gauss-seidel-coefficient-form}). For convenience I have repeated it here below:

$$
x_i^{n+1} = \frac{1}{a_{ii}}\left(b_i - \sum_{j=1}^{i-1}l_{ij}x_j^{n+1} - \sum_{j=i+1}^{N}u_{ij}x_j^{n}\right)
$$

The right-hand side is now updating the left-hand side fully, but we want to write the right-hand side in the form of Eq.(\ref{eq:under-relaxation}). We can do this by introducing the relaxation parameter $\omega$ (as it is usually used in the literature) and write:

$$
x_i^{n+1} = \omega\frac{1}{a_{ii}}\left(b_i - \sum_{j=1}^{i-1}l_{ij}x_j^{n+1} - \sum_{j=i+1}^{N}u_{ij}x_j^{n}\right) + (1-\omega ) x_i^{n}
$$

Or, rearranging terms, we can also write:

$$
x_i^{n+1} = (1-\omega ) x_i^{n} + \frac{\omega}{a_{ii}}\left(b_i - \sum_{j=1}^{i-1}l_{ij}x_j^{n+1} - \sum_{j=i+1}^{N}u_{ij}x_j^{n}\right)
\tag{eq:gauss-seidel-sor-coefficient-form}
$$

This is perhaps more commonly used in the literature, but both are the same equation. So, we multiply the term that updates $\mathbf{x}^{n+1}$ by $\omega$ and then, we have to introduce $\mathbf{x}^n$ and multiply that by $(1-\omega)$ to get a weighted average between the solution at $n$ and $n+1$. The question is, what values can we choose for $\omega$ so that the solution is still stable?

Well, we have already established that for the simulation to converge, we require that the spectral radius of the matrix $\mathbf{G}$ is less than 1. As long as that is given, we will get convergence. The spectral radius, as defined in Eq.(\ref{eq:spectral-radius}), says that it is the largest (absolute) eigenvalue of a matrix. So, we need to do 2 things: First, we need to derive the matrix $\mathbf{G}$ which contains $\omega$, and then we need to compute its eigenvalues, and solve those for $\omega$.

Let's start by deriving the iteration matrix $\mathbf{G}$. From Eq.(\ref{eq:comparison-jacobi-gs-2}), we saw that the Gauss-Seidel method is given as:

$$
(\mathbf{D} + \mathbf{L})\mathbf{x}^{n+1} + \mathbf{U}\mathbf{x}^{n}=\mathbf{b}
$$

We bring the contributions of $n$ to the right-hand side and get:

$$
(\mathbf{D} + \mathbf{L})\mathbf{x}_{GS}^{n+1} = \mathbf{b} - \mathbf{U}\mathbf{x}^{n}
\tag{eq:gauss-seidel-solved-n1}
$$

Here, I have started to use the notation $\mathbf{x}_{GS}^{n+1}$ to show that this unknown solution vector is obtained from the Gauss-Seidel method. Now, let's remind ourselves of the SOR definition. We said that this is coming from Eq.(\ref{eq:under-relaxation}), and so we can write this, using the SOR factor $\omega$, instead of $\alpha$, as:

$$
\mathbf{x}^{n+1}=\omega\mathbf{x}_{GS}^{n+1} + (1-\omega)\mathbf{x}^n 
\tag{eq:sor-definition}
$$

Here, $\mathbf{x}_{GS}^{n+1}$ is the solution we obtain from the Gauss-Seidel algorithm, i.e. Eq.(\ref{eq:gauss-seidel-solved-n1}). Let's isolate that on the left-hand side:

$$
\omega\mathbf{x}_{GS}^{n+1} = \mathbf{x}^{n+1} - (1-\omega)\mathbf{x}^n
$$

We divide this now by $\omega$, which provides us with the final form:

$$
\mathbf{x}_{GS}^{n+1} = \frac{1}{\omega}\left(\mathbf{x}^{n+1} - (1-\omega)\mathbf{x}^n\right)
$$

We can insert this definition now back into Eq.(\ref{eq:gauss-seidel-solved-n1}) and get:

$$
(\mathbf{D} + \mathbf{L})\frac{1}{\omega}\left(\mathbf{x}^{n+1} - (1-\omega)\mathbf{x}^n\right) = \mathbf{b} - \mathbf{U}\mathbf{x}^{n}
$$

We multiply this equation by $\omega$ and get:

$$
(\mathbf{D} + \mathbf{L})\left(\mathbf{x}^{n+1} - (1-\omega)\mathbf{x}^n\right) = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n}
$$

Let us now expand the left-hand side. We have:

$$
(\mathbf{D} + \mathbf{L})\mathbf{x}^{n+1} - (\mathbf{D} + \mathbf{L})(1-\omega)\mathbf{x}^n = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n}
$$

Let's further expand the left-hand side:

$$
\mathbf{D}\mathbf{x}^{n+1} + \mathbf{L}\mathbf{x}^{n+1} - \mathbf{D}(1-\omega)\mathbf{x}^n - \mathbf{L}(1-\omega)\mathbf{x}^n = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n}
$$

Let me rewrite this equation as:

$$
\mathbf{D}\mathbf{x}^{n+1} - \mathbf{D}(1-\omega)\mathbf{x}^n + \mathbf{L}\mathbf{x}^{n+1} - \mathbf{L}(1-\omega)\mathbf{x}^n = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n}
\tag{eq:gs-derivation-1}
$$

Let's look at Eq.(\ref{eq:sor-definition}) for a moment. If we multiply this equation on each side with $\mathbf{L}$, then we get:

$$
\mathbf{L}\mathbf{x}^{n+1}=\mathbf{L}\omega\mathbf{x}_{GS}^{n+1} + \mathbf{L}(1-\omega)\mathbf{x}^n
$$

Solving for $\mathbf{x}_{GS}^{n+1}$, we get:

$$
\mathbf{L}\omega\mathbf{x}_{GS}^{n+1} = \mathbf{L}\mathbf{x}^{n+1} - \mathbf{L}(1-\omega)\mathbf{x}^n 
$$

Noting that $\mathbf{x}_{GS}^{n+1}=\mathbf{x}^{n+1}$, since we already inserted this into our Gauss-Seidel algorithm, we can re-write Eq.(\ref{eq:gs-derivation-1}) as:

$$
\mathbf{D}\mathbf{x}^{n+1} - \mathbf{D}(1-\omega)\mathbf{x}^n + \mathbf{L}\omega\mathbf{x}^{n+1} = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n}
\tag{eq:gs-derivation-2}
$$

Now, you may be saying, we can do the same for the first two terms involving $\mathbf{D}$, right? OK, let's do that and see what happens. We can write the SOR definition, now multiplied by $\mathbf{D}$, as:

$$
\mathbf{D}\omega\mathbf{x}_{GS}^{n+1} = \mathbf{D}\mathbf{x}^{n+1} - \mathbf{D}(1-\omega)\mathbf{x}^n 
$$

We insert that definition into Eq.(\ref{eq:gs-derivation-2}), noting again that $\mathbf{x}_{GS}^{n+1}=\mathbf{x}^{n+1}$, and we get:

$$
\mathbf{D}\omega\mathbf{x}^{n+1} + \mathbf{L}\omega\mathbf{x}^{n+1} = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n}
\tag{eq:gs-derivation-3}
$$

Now we notice that each term contains $\omega$, so we can divide each term by it:

$$
\mathbf{D}\mathbf{x}^{n+1} + \mathbf{L}\mathbf{x}^{n+1} = \mathbf{b} - \mathbf{U}\mathbf{x}^{n}
$$

Finally, we can simplify the left-hand side to:

$$
(\mathbf{D} + \mathbf{L})\mathbf{x}^{n+1} = \mathbf{b} - \mathbf{U}\mathbf{x}^{n}
$$

We just obtained Eq.(\ref{eq:gauss-seidel-solved-n1}) again, which was the starting point for our derivation. So, surely, this is not what we want. Instead, let's go back to Eq.(\ref{eq:gs-derivation-2}) and leave the terms involving $\mathbf{D}$ as they are on the left-hand side. The next step in our derivation requires us to bring $\mathbf{D}(1-\omega)\mathbf{x}^n$ to the right-hand side. This gives us:

$$
\mathbf{D}\mathbf{x}^{n+1} + \mathbf{L}\omega\mathbf{x}^{n+1} = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n} + \mathbf{D}(1-\omega)\mathbf{x}^n
$$

Now we introduce a little trick to make the signs the same. The following definition holds:

$$
(1-\omega)=-(\omega - 1)
$$

If you don't believe me, insert any value for $\omega$ and you will see this is correct. We can use this to rewrite our equation as:

$$
\mathbf{D}\mathbf{x}^{n+1} + \mathbf{L}\omega\mathbf{x}^{n+1} = \omega\mathbf{b} - \omega\mathbf{U}\mathbf{x}^{n} - \mathbf{D}(\omega - 1)\mathbf{x}^n
$$

The second and third term both contain a $\mathbf{x}^n$ on the right-hand side, let's factor that out and write:

$$
\mathbf{D}\mathbf{x}^{n+1} + \mathbf{L}\omega\mathbf{x}^{n+1} = \omega\mathbf{b} - (\omega\mathbf{U} - \mathbf{D}(\omega - 1))\mathbf{x}^n
$$

Now, onto the left-hand side, we can rewrite this in a more compact form as:

$$
(\mathbf{D} + \omega\mathbf{L})\mathbf{x}^{n+1} = \omega\mathbf{b} - (\omega\mathbf{U} - \mathbf{D}(\omega - 1))\mathbf{x}^n
\tag{eq:gs-sor-final}
$$

Comparing this with the original Gauss-Seidel method without any SOR treatment, we had:

$$
(\mathbf{D} + \mathbf{L})\mathbf{x}_{GS}^{n+1} = \mathbf{b} - \mathbf{U}\mathbf{x}^{n}
$$

We see that the left-hand side remains the same, apart from the relaxation factor $\omega$. We also see that the right-hand side is slightly modified, but more crucially, we see that dividing this entire equation by $\omega$ will not get rid of it. This is, perhaps, the best intuitive approach I can show you that this form is correct, and not the one we saw before where $\omega$ cancels from the equation. 

We can now, again, define an exact solution $\mathbf{x}^*$ and form the same SOR system:

$$
(\mathbf{D} + \omega\mathbf{L})\mathbf{x}^* = \omega\mathbf{b} - (\omega\mathbf{U} - \mathbf{D}(\omega - 1))\mathbf{x}^*
\tag{eq:gs-sor-exact}
$$

If we subtract Eq.(\ref{eq:gs-sor-exact}) from Eq.(\ref{eq:gs-sor-final}), we get:

$$
(\mathbf{D} + \omega\mathbf{L})(\mathbf{x}_{GS}^{n+1} - \mathbf{x}^*) = (\omega\mathbf{U} - \mathbf{D}(\omega - 1))(\mathbf{x}^{n} - \mathbf{x}^*)
$$

Introducing the error notation again where $\mathbf{x}_{GS}^{n+1} - \mathbf{x}^*=e^{n+1}$ and $\mathbf{x}^{n} - \mathbf{x}^* = e^n$, we get:

$$
(\mathbf{D} + \omega\mathbf{L})e^{n+1} = (\omega\mathbf{U} - \mathbf{D}(\omega - 1))e^n
$$

Solving for the error $e^{n+1}$, we get:

$$
e^{n+1} = (\mathbf{D} + \omega\mathbf{L})^{-1}(\omega\mathbf{U} - \mathbf{D}(\omega - 1))e^n
$$

Thus, we can define our matrix $\mathbf{G}$ as:

$$
\mathbf{G} = (\mathbf{D} + \omega\mathbf{L})^{-1}(\omega\mathbf{U} - \mathbf{D}(\omega - 1))
$$

From our spectral radius definition, we know that this converges for $\rho(\mathbf{G})\lt 1$. However, we also know that the spectral radius is nothing more than the largest (absolute) eigenvalue of $\mathbf{G}$, i.e. $\rho(\mathbf{G})=|\lambda_i|$. So, we need the eigenvalues of $\mathbf{G}$. The eigenvalues are obtained by evaluating $\text{det}(\mathbf{G}-\lambda\mathbf{I})$ and setting this to zero, i.e.:

$$
\text{det}(\mathbf{G}-\lambda\mathbf{I})=0
$$

Let's apply that to an example for the matrix:

$$
\mathbf{G}=
\begin{bmatrix}
2 & 1\\
1 & 2
\end{bmatrix}
$$

The eigenvalues can be computed as:

$$
\text{det}(\mathbf{G}-\lambda\mathbf{I})=
\begin{bmatrix}
2-\lambda & 1\\
1 & 2-\lambda
\end{bmatrix}=
\lambda^2 - 4\lambda + 3 = a\lambda^2 + b\lambda + c = 0
$$

Here, I have defined the constants $a=1$, $b=-4$, and $c=3$. We can use the abc formula to solve this quadratic equation as:

$$
\lambda_{1,2}=\frac{-b\pm\sqrt{b^2-4ac}}{2a} = \frac{-(-4)\pm\sqrt{(-4)^2 - 4(1)(3)}}{2(1)}=\frac{4\pm\sqrt{16 - 12}}{2}=\frac{4\pm 2}{2}\\[1em]
\lambda_1 = \frac{4 + 2}{2} = \frac{6}{2} = 3 \\[1em]
\lambda_2 = \frac{4 - 2}{2} = \frac{2}{2} = 1
$$

Since $\text{det}(\mathbf{G}-\lambda\mathbf{I})$ is set equal to zero when evaluating the eigenvalues, we also have $\lambda^2 - 4\lambda + 3 = 0$ as we saw above. We can rewrite this equation, knowing the eigenvalues now, as:

$$
\lambda^2 - 4\lambda + 3 = (3-\lambda)(1-\lambda)
$$

If we insert either 1 or 3 into the equation, i.e. our eigenvalues, then the right-hand side will be equal to zero. We can generalise this as:

$$
\lambda^2 - 4\lambda + 3 = (\lambda_1 - \lambda)(\lambda_2 - \lambda)
$$

More generally, we can write:

$$
\text{det}(\mathbf{G}-\lambda\mathbf{I}) = (\lambda_1 - \lambda)(\lambda_2 - \lambda) ... (\lambda_n - \lambda)
$$

Remember what our goal is. We have the matrix $\mathbf{G}$ given, and we want to relate that to our eigenvalues $\lambda_1,\lambda_2,...,\lambda_n$, specifically, we want to know which is the largest eigenvalue, so that we can evaluate the spectral radius. In the equation given above, $\lambda$ is just a variable, whereas $\lambda_i$ are our eigenvalues. Therefore, we can set $\lambda=0$ and get:

$$
\text{det}(\mathbf{G}-0\mathbf{I}) = (\lambda_1 - 0)(\lambda_2 - 0) ... (\lambda_n - 0)\\[1em]
\text{det}(\mathbf{G}) = \lambda_1\lambda_2 ... \lambda_n\\[1em]
\text{det}(\mathbf{G}) = \prod_{i=1}^n\lambda_i
\tag{eq:characteristic-polynomial}
$$

Here, $n$ is the number of rows and columns of $\mathbf{G}$, i.e. we have as many eigenvalues as there are rows in the matrix. How does this help us? Well, we said that in order to have convergence, we need the spectral radius, that is, the largest eigenvalue, to be smaller than 1. However, Eq.(\ref{eq:characteristic-polynomial}) gives the product of all eigenvalues, we can't easily extract the largest eigenvalue from it.

Here is the trick: We say that in order for the sepctral radius to be larger than 1, we need to have at least one eigenvalue that is larger than 1. If all eigenvalues are less than 1, then their product will also be less then one. However, if one eigenvalue is larger than 1, it is possible for the product to be larger than 1 as well.

Now, I should state here that this is only one condition, for it is possible to have eigenvalues that are larger than 1, but the products of all eigenvalues is less then one. If I have $\lambda_1=5$ and $\lambda_2=0.1$, then the product is $\lambda_1\lambda_2=5\cdot 0.1 = 0.5$. The product is less then one, but the largest individual eigenvalue is larger than one.

So, while the Eq.(\ref{eq:characteristic-polynomial}) doesn't say anything about the largest eigenvalues, we can still use it. That is, if the product of all eigenvalues is larger than 1, that tells us that at least one eigenvalue is larger than one as well, as already established. It turns out that, even though this is a weaker requirement than checking the absolute largest eigenvalue, it still works in practice very well.

So, what we do is to say that in order for the SOR Gauss-Seidel method to work and coonverge, we have the requirement:

$$
\text{det}(\mathbf{G}) = \prod_{i=1}^n|\lambda_i| \lt 1
\tag{eq:det-lambda-equivalence}
$$

Notice that it is the absolute value of the eigenvalues now, i.e. the eigenvalues can be negative. This comes from the definition of the spectral radius, i.e. Eq.(\ref{eq:spectral-radius}).

OK, so let's evaluate this expression. First, let us remind ourselves what the matrix $\mathbf{G}$ was, which I have repeated below for convenience:

$$
\mathbf{G} = (\mathbf{D} + \omega\mathbf{L})^{-1}(\omega\mathbf{U} - \mathbf{D}(\omega - 1))
$$

We want to take the determinant of this now. There are a few rules for determinants that we can apply here. The following [properties](https://en.wikipedia.org/wiki/Determinant#Multiplicativity_and_matrix_groups) hold:

$$
\text{det}(\mathbf{AB})=\text{det}(\mathbf{A})\text{det}(\mathbf{B})\\[1em]
\text{det}(\mathbf{A}^{-1})=\frac{1}{\text{det}(\mathbf{A})}
$$

Using these properties, we can write the determinant of $\mathbf{G}$ as:

$$
\text{det}(\mathbf{G}) = \text{det}\left[(\mathbf{D} + \omega\mathbf{L})^{-1}(\omega\mathbf{U} - \mathbf{D}(\omega - 1))\right]\\[1em]
\text{det}(\mathbf{G}) = \text{det}\left[(\mathbf{D} + \omega\mathbf{L})^{-1}\right]\text{det}\left[\omega\mathbf{U} - \mathbf{D}(\omega - 1)\right]\\[1em]
\text{det}(\mathbf{G}) = \frac{\text{det}\left[\omega\mathbf{U} - \mathbf{D}(\omega - 1)\right]}{\text{det}\left[\mathbf{D} + \omega\mathbf{L}\right]}
\tag{eq:determinant-of-G}
$$

We now need to evaluate these determinants separately. For a 3-by-3 matrix, we can use the rule of Sarrus, where we copy the first two columns to the right, and develop our determinant as shown in the figure below on the left:

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\determinant_decomposition.png" alt="Development of the determinant for a dense, diagonal, and lower, upper triangualr matrix, where the lower and upper matrices have a determinant of zero." class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

Here, we have that the determinant for the matrix $\mathbf{A}$ is:

$$
\text{det}(\mathbf{A})=aei - gec + bfg - hfa + cdh - idb
$$

We can see that the determinant essentially adds and subtracts products of diagonals and off-diagonals. Let's apply this to the diagonal matrix $\mathbf{D}$. This matrix has, by definition, the main diagonal of matrix $\mathbf{A}$. As long as the values on the diagonal are not zero (which they are not for CFD applications), then we will have a non-zero determinant. This is also shown in the figure above.

In the same figure, I have also given the determinants for the lower and upper triangular matrix $\mathbf{L}$ and $\mathbf{U}$. I have removed entries here which would be zero. We can see that if we apply the rule of Sarrus to either $\mathbf{L}$ or $\mathbf{U}$, there will always be at least one term that will be zero. As a result, both $\mathbf{L}$ and $\mathbf{U}$ will have a determinant that is zero.

We can use that now in Eq.(\ref{eq:determinant-of-G}). Let us look at the determinants for both the numerator and denominator separately. For the numerator, we have:

$$
\text{det}\left[\omega\mathbf{U} - \mathbf{D}(\omega - 1)\right]
$$

If we look now at diagonal entries, only $\mathbf{D}$ has one full diagonal with non-zero entries. Therefore, we can write the determinant expression as:

$$
\text{det}\left[\omega\mathbf{U} - \mathbf{D}(\omega - 1)\right]=\sum_{i=1}^{n}d_{ii}(\omega - 1)
$$

Here, we are summing the individual contributions in $\mathbf{D}$, i.e. $d_{ii}$ with the expression given inside the brackets, i.e. $(\omega - 1)$, which gives us the only non-zero diagonal contribution and thus determinant. We can remove the summation and write this in matrix form as:

$$
\sum_{i=1}^{n}d_{ii}(\omega - 1) = \text{det}(\mathbf{D})(\omega - 1)^n
$$

Since the determinant of $\mathbf{D}$ is just the sum of all of its diagonal entries (as all other off-diagonals will be zero), we can replace the summation here by the determinant. The only thing we have to change here is that the expression given in parenthesis is now raised to the power of $n$, as we have to multiply it $n$-times with the coefficients in $\mathbf{D}$.

OK, return to Eq.(\ref{eq:determinant-of-G}), the denominator can be evaluated as:

$$
\text{det}\left[\mathbf{D} + \omega\mathbf{L}\right] = \text{det}(\mathbf{D})
$$

Here, we use the same argument again that only $\mathbf{D}$ will contribute to non-zero products in the development of the diagonal, whereas all terms of the lower triangular matrix will have zeros in their product and thus it plays no further role ion the development of the determinant.

If we insert now the simplified expressions for the numerator and denominator back into Eq.(\ref{eq:determinant-of-G}), we get:

$$
\text{det}(\mathbf{G}) = \frac{\text{det}\left[\omega\mathbf{U} - \mathbf{D}(\omega - 1)\right]}{\text{det}\left[\mathbf{D} + \omega\mathbf{L}\right]}\\[1em]
\text{det}(\mathbf{G}) = \frac{\text{det}(\mathbf{D})(\omega - 1)^n}{\text{det}(\mathbf{D})}\\[1em]
\text{det}(\mathbf{G}) = \frac{(\omega - 1)^n}{1}\\[1em]
\text{det}(\mathbf{G}) = (\omega - 1)^n
\tag{eq:determinant-of-G-2}
$$

We can use this expression now and insert that back into Eq.(\ref{eq:det-lambda-equivalence}). This gives us:

$$
(\omega - 1)^n = \prod_{i=1}^n|\lambda_i| \lt 1\\[1em]
|(\omega - 1)^n| \lt 1
$$

The determinant can be positive or negative, but since it is used here to compute eigenvalues (for which we have to take the absolute value), we are now using the absolute value of $(\omega - 1)^n$. We need to solve this expression now for $\omega$. The first step is to take the $n^{th}$-root ($\sqrt[n]{}$) of both side. This gives us:

$$
|(\omega - 1)| \lt 1
$$

Now we can check for which ranges of $\omega$ the inequality holds. We can that if $\omega=0$, the expression becomes:

$$
|(0 - 1)| \lt 1\\
|-1| \lt 1\\
1 \lt 1
$$

Thus, $\omega=0$ is the lower limit at which point convergence is no longer given (it is also the trivial answer as an under-relaxation factor of zero just means we are not updating the solution at all!). So, we need to set $\omega\gt 0$ to get stable convergence. If we set $\omega=2$, then we get:

$$
|(2 - 1)| \lt 1\\
|1| \lt 1\\
1 \lt 1
$$

Again, a value of $\omega=2$ is just over the limit and we will get neither convergence nor divergence. So, the upper limit of $\omega$ is two, and we can write this in a compact form as:

$$
0\lt \omega \lt 2
$$

This is the requirement for $\omega$, and, unlike under-relaxation factors, in the SOR method (which literally stands for successive over-relaxation), we are able to over-relax the Gauss-Seidel method and thus accelerate convergence. This is one rare instance in CFD where you get quite a lot of return on essentiually a zero investment (all we have to do is to multiply our Gauss-Seidel algorithm with some $\omega$s).

Thus, in the extreme case of $\omega \approx 2$, we can get a convergence acceleration of 4 compared to the Jacobi method. Not bad for simply setting ```xOld=xNew``` and multiplying by one factor.

Damn, I just realised, I'll be flying out of heathrow in a few weeks again. Oh boy ...

### Lower-Upper Symmetric Gauss-Seidel (LU-SGS)

If you thought we would be done with Gauss, I'm sorry to dissapoint you, he is here to stay and we need to talk about yet another method. It is a bit of a stretch to call it *another method*. The lower-upper, syemmtric Gauss-Seidel, or LU-SGS method, is nothing else than just the pure symmetric Gauss-Seidel method. It seems perhaps more of a footnote lost to history that some people prefer to use the terminology LU-SGS, while others simply stick with the symmetric Gauss-Seidel terminology.

To the best of my knowledge, LU-SGS is what is what people use in the compressible flow community, while others (especially incompressible flow developers) use the symemtric Gauss-Seidel terminology. Again, both methods are algebraically identical, they just have a different name.

However, despite being algebraically different, the way the methods are derived are slightly different. For the Gauss-Seidel method, we simply set $x_i^n$ to $x_i^{n+1}$ if the solution was already updated inside the loop from $i=1$ to $N$. That's not much of a derivation, granted.

In the LU-SGS method, the starting point is, in a very anticlimatic fashion, the LU decomposition. So, let us look at the derivation and see for ourselves that the LU-SGS method does indeed simply produce the symmetric Gauss-Seidel method.

The starting point for the LU decomposition was to split the matrix $\mathbf{A}$ into a lower and upper triangular matrix, such that if we multiplied both, we would get the original matrix $\mathbf{A}$ back. We can express this as:
 
$$
\mathbf{A} = \tilde{\mathbf{L}} \, \tilde{\mathbf{U}}
$$

We saw in the section on the LU decomposition that forming $\tilde{\mathbf{L}}$ and $\tilde{\mathbf{U}}$ costs just as much as doing the Gauss elimination, which is prohibitively slow. For that reason, we don't want to find an exact LU decomposition, but rather an approximate LU decomposition. Ideally, we only want to use information we already have as well, so we don't need to compute extra bits and pieces.

Thus far in our discussion, we have often decomposed the coefficient matric $\mathbf{A}$ as:

$$
\mathbf{A} = \mathbf{D} + \mathbf{L} + \mathbf{U}
$$

Since we have all matrices available, our goal is to use these to find an approximate LU decomposition. A first idea may be to use:

$$
(\mathbf{D} + \mathbf{L})(\mathbf{D} + \mathbf{U})
$$

However, if we carry out the multiplication between the brackets, then we obtained the following:

$$
(\mathbf{D} + \mathbf{L})(\mathbf{D} + \mathbf{U})=\mathbf{D}^2+\mathbf{D}\mathbf{U}+\mathbf{L}\mathbf{D}+\mathbf{L}\mathbf{U}
$$

We now get a term $\mathbf{D}^2$, as well as a combined $\mathbf{L}\mathbf{U}$. Why is this a problem. Remember, the goal is to express $\tilde{\mathbf{L}} \, \tilde{\mathbf{U}}$ through already know matrices, and we said that when we decompose the coefficient matrix $\mathbf{A}$, we did so by writing $\mathbf{D} + \mathbf{L} + \mathbf{U}$. So, we are looking for approximation to $\tilde{\mathbf{L}} \, \tilde{\mathbf{U}}$ that result in $\mathbf{D} + \mathbf{L} + \mathbf{U}$.

Since $(\mathbf{D} + \mathbf{L})(\mathbf{D} + \mathbf{U})$ by itself gives us additional terms, such as $\mathbf{D}^2$, we need to modify this term slightly. Now think about it yourself. If you have the number $4^2$, and you want to get rid of the power, you could take the square root, yes, but that wouldn't work on matrices. But, we can write $4^2=4\cdot 4$, and so, if we now want to get rid of the power of two, we could also write $(4\cdot 4)/4$ or $4^2\cdot 4^{-1}=4$.

We can take the inverse of a matrix, and so, if we want to get rid of power in the $\mathbf{D}^2$, because we are trying to obtain a form of $\mathbf{D} + \mathbf{L} + \mathbf{U}$, then we simply multiply by $\mathbf{D}^{-1}$, because $\mathbf{D}^2\mathbf{D}^{-1}=\mathbf{D}$. Makes sense? So, let's do that. We multiply both parenthesis by $\mathbf{D}^{-1}$ and we get:

$$
(\mathbf{D} + \mathbf{L}) \mathbf{D}^{-1} (\mathbf{D} + \mathbf{U})
\tag{eq:lusgs-1}
$$

Let's carry out the multiplication for the first term with $\mathbf{D}^{-1}$. We get:

$$
(\mathbf{D} + \mathbf{L}) \mathbf{D}^{-1}=\mathbf{D}\mathbf{D}^{-1} + \mathbf{L}\mathbf{D}^{-1} = \mathbf{I} + \mathbf{L}\mathbf{D}^{-1}
\tag{eq:lusgs-2}
$$

If our diagonal contains entries like $d_{ii}$, then its inverse ($\mathbf{D}^{-1}$) will contain entries like $1/d_{ii}$. Remember, if you only have the diagonal in a matrix (and, after all, $\mathbf{D}$ is still a matrix, but the letter D indicates that it only has entries on thje diagonal), we can find the inverse to that matrix with ease by inverting the diagonal.

Thus, if we compute now $\mathbf{D}\mathbf{D}^{-1}$, then we will evaluate products like $d_{ii}(1/d_{ii})=d_{ii}/d_{ii}=1$, and therefore, $\mathbf{D}\mathbf{D}^{-1}=\mathbf{I}$, where $\mathbf{I}$ is the identity matrix.

We can now multiply the term obtained from Eq.(\ref{eq:lusgs-2}) with the second term in Eq.(\ref{eq:lusgs-1}) and get:

$$
(\mathbf{I} + \mathbf{L}\mathbf{D}^{-1})(\mathbf{D} + \mathbf{U}) = \mathbf{I}\mathbf{D} + \mathbf{I}\mathbf{U} + \mathbf{L}\mathbf{D}^{-1}\mathbf{D} + \mathbf{L}\mathbf{D}^{-1}\mathbf{U}
$$

Acknowledging that any matrix multiplied by the identity matrix will give back the original matrix. Furthermore, we see that we have $\mathbf{D}^{-1}\mathbf{D}$ again, which is the identity matrix. Since the identity matrix does not modify a matrix if we multiply by it, we can drop it can obtain:

$$
\underbrace{\mathbf{D} + \mathbf{U} + \mathbf{L}}_\mathbf{A} + \mathbf{L}\mathbf{D}^{-1}\mathbf{U}
$$

We can write this as:

$$
\mathbf{A} - \mathbf{L}\mathbf{D}^{-1}\mathbf{U} = \mathbf{D} + \mathbf{L} + \mathbf{U}
$$

We have almost reached our goal, that is, to find an approximation to $\tilde{\mathbf{L}} \, \tilde{\mathbf{U}}$ that results in a decomposition of the form $\mathbf{D} + \mathbf{L} + \mathbf{U}$, but we have obtained an additional $\mathbf{L}\mathbf{D}^{-1}\mathbf{U}$ term, which we now need to get rid off. How do we do that?

Well, we could turn to our friends at Heathrow airport and asked a border patrol agent what they would do, but their answer probably would be just *deport it*. Now, I have checked, but there doesn't seem to be a *deportation* operator we could use here, and so we have to use a bit of logic instead. First, we write the term on the right-hand side so that the matrix is again isolated on the left. We have:

$$
\mathbf{A} = (\mathbf{D} + \mathbf{L} + \mathbf{U}) + \mathbf{L}\mathbf{D}^{-1}\mathbf{U}
\tag{eq:lusgs-with-additional-term}
$$

Now comes our logic, or trick. We say: $\mathbf{L}\mathbf{D}^{-1}\mathbf{U}\ll \mathbf{D} + \mathbf{L} + \mathbf{U}$ and therefore we drop it altogether. Why can we do this? Well, we can use the so-called *boundness* property, which was formulated in 1958 by Scarborough. To do that, we need to write $\mathbf{L}\mathbf{D}^{-1}\mathbf{U}$ in coefficient form first. For that, we make the following substitutions: 

$$
\mathbf{D}^{-1}=\frac{1}{d_{ii}}\\[1em] 
\mathbf{L}=\sum_{j=1}^{i-1}l_{ij}\\[1em] 
\mathbf{U}=\sum_{j=i+1}^{N}u_{ij}\\[1em]
$$

Then, we can write:

$$
\mathbf{L}\mathbf{D}^{-1}\mathbf{U} = \frac{\left(\sum_{j=1}^{i-1}l_{ij} \right)\left(\sum_{j=i+1}^{N}u_{ij} \right)}{d_{ii}}
\tag{eq:LDU}
$$

Scaraborough says that the following inequality must hold:

$$
\frac{\sum_{j=1}^{i-1}l_{ij} + \sum_{j=i+1}^{N}u_{ij}}{d_{ii}}=
\begin{cases}
\le 1 & \text{for all cells}\\[1em]
\lt 1 & \text{in at least one cell}
\end{cases}
$$

Compare the two states. While the boundness property adds the summations in the numerator, we multiply them in Eq.(\ref{eq:LDU}). The boundness property is a necessary condition for condition for convergence, and, indeed, the diagonal typically contains the largest elements in the entire coefficient matrix. Thus, if we look at $\mathbf{L}\mathbf{D}^{-1}\mathbf{U}$, we see that we scale both $\mathbf{L}$ and $\mathbf{U}$ by the inverse of $\mathbf{D}$, and $\mathbf{D}$ must be greater than $\mathbf{L}$ and $\mathbf{U}$ combined.

Why is this typically the case? The timestep only appears in the $\mathbf{D}$ matrix, not in the $\mathbf{L}$ and $\mathbf{U}$ matrix, and this typically results $\mathbf{D}$ being at least an order of magnitude greater. For example, if we had the following PDE:

$$
\frac{\partial T}{\partial t} = \alpha\frac{\partial^2 T}{\partial x^2}
$$

Then we could discretise it with a simple first-order forward in time and second-order central scheme in space (using the finite difference method) as:

$$
\frac{T^{n+1}_i-T^n_i}{\Delta t} = \alpha\frac{T^{n+1}_{i+1} - 2T^{n+1}_i + T^{n+1}_{i-1}}{(\Delta x)^2}
$$

We can now write this in terms of the factors of $T$ and get:

$$
T^{n+1}_{i}\left[\frac{1}{\Delta t}\right] + T^{n}_{i}\left[\frac{-1}{\Delta t}\right] = T^{n+1}_{i+1}\left[\frac{\alpha}{(\Delta x)^2}\right] + T^{n+1}_{i}\left[\frac{-2\alpha}{(\Delta x)^2}\right] + T^{n+1}_{i-1}\left[\frac{\alpha}{(\Delta x)^2}\right]
$$

Now we collect terms at $n+1$ on the left hand side and terms of $n$ on the right-hand side:

$$
T^{n+1}_{i}\left[\frac{1}{\Delta t}\right] + T^{n+1}_{i+1}\left[\frac{-\alpha}{(\Delta x)^2}\right] + T^{n+1}_{i}\left[\frac{2\alpha}{(\Delta x)^2}\right] + T^{n+1}_{i-1}\left[\frac{-\alpha}{(\Delta x)^2}\right]= T^{n}_{i}\left[\frac{1}{\Delta t}\right]
$$

We can now combine the two terms at $T^{n+1}_{i}$:

$$
T^{n+1}_{i}\left[\frac{1}{\Delta t} + \frac{2\alpha}{(\Delta x)^2}\right] + T^{n+1}_{i+1}\left[\frac{-\alpha}{(\Delta x)^2}\right] + T^{n+1}_{i-1}\left[\frac{-\alpha}{(\Delta x)^2}\right]= T^{n}_{i}\left[\frac{1}{\Delta t}\right]
$$

Terms at $i$ will go into our diagonal matrix $\mathbf{D}$, while terms at $i+1$ and $i-1$ will go into $\mathbf{U}$ and $\mathbf{L}$, respectively. Thus, we could write this as:

$$
T^{n+1}_{i}\mathbf{D} + T^{n+1}_{i+1}\mathbf{U} + T^{n+1}_{i-1}\mathbf{L}= T^{n}_{i}\left[\frac{1}{\Delta t}\right]
$$

Therefore, we have:

$$
\mathbf{D} = \left[\frac{1}{\Delta t} + \frac{2\alpha}{(\Delta x)^2}\right]\\[1em]
\mathbf{U} = \left[\frac{-\alpha}{(\Delta x)^2}\right]\\[1em]
\mathbf{L} = \left[\frac{-\alpha}{(\Delta x)^2}\right]
\tag{eq:LDU2}
$$

Now we can do some order of magnitude estimates. Let's start with $\alpha$. Looking at some values for the [thermal diffusivity](https://en.wikipedia.org/wiki/Thermal_diffusivity#Thermal_diffusivity_of_selected_materials_and_substances), we have values between $10^{-8}\lt \alpha \lt 10^{-4}$. For example, we have $\alpha = 1.9\cdot 10^{-5}$ for air and $\alpha = 1.47 \cdot 10^{-7}$ for water (measured in $m^2/s$).

From the CFL condition for purely diffusion flows, we have:

$$
CFL = \frac{\alpha\Delta t}{(\Delta x)^2}
$$

Let's rewrite this as:

$$
\Delta t = \frac{CFL}{\alpha}(\Delta x)^2
$$

If we say that $\alpha$ ranges between $10^{-8}\lt \alpha \lt 10^{-4}$, let's pick an average $\alpha=10^{-6}$. Let's insert this into our equation:

$$
\Delta t = \frac{CFL}{10^{-6}}(\Delta x)^2\\[1em]
\Delta t = 10^6 CFL(\Delta x)^2
$$

We can now assume that our CFL number is typically of the order of 1 to 100. Then, we can write:

$$
10^6 (\Delta x)^2 \le \Delta t \le 10^8 (\Delta x)^2
$$

In other words, it is fair to say that $\Delta t \gg \Delta x$, and this is true for convection-driven flows as well (in this case we only looked at thermal diffusion). Let's pick a moderate CFL number of 10. Let's also pick an arbitrary spacing of $\Delta x = 10^{-5}$. Then, we can compute the timestep to be:

$$
\Delta t = \frac{CFL}{10^{-6}}(\Delta x)^2 = \frac{10}{10^{-6}}(10^{-5})^2 = 10\cdot 10^{6} \cdot 10^{-10} = 10^{-3}
$$

Thus, we have $\Delta t=10^{-3}$ while we have $\Delta x = 10^{-5}$, and indeed, we have $\Delta t \gg \Delta x$. Let's go back to Eq.(\ref{eq:LDU2}) and insert these values:

$$
\mathbf{D} = \left[\frac{1}{\Delta t} + \frac{2\alpha}{(\Delta x)^2}\right] = \left[\frac{1}{10^{-3}} + \frac{2\cdot 10^{-6}}{10^{-10}}\right] = \left[1000 + 20000\right] = 21000 \\[1em]
\mathbf{U} = \left[\frac{-\alpha}{(\Delta x)^2}\right] = \left[\frac{-10^{-6}}{(\Delta 10^{-5})^2}\right] = \left[\frac{-10^{-6}}{10^{-10}}\right] = -10000 \\[1em]
\mathbf{L} = \left[\frac{-\alpha}{(\Delta x)^2}\right] = \left[\frac{-10^{-6}}{(\Delta 10^{-5})^2}\right] = \left[\frac{-10^{-6}}{10^{-10}}\right] = -10000
$$

Let's go back to our original statement that $\mathbf{L}\mathbf{D}^{-1}\mathbf{U}\ll \mathbf{D} + \mathbf{L} + \mathbf{U}$. We can now write out our order of magnitude estimates for this particular case:

$$
\mathbf{L}\mathbf{D}^{-1}\mathbf{U}\ll \mathbf{D} + \mathbf{L} + \mathbf{U}\\
|-10000||-10000|/21000 \ll 21000 + |-10000| + |-10000|\\
4762 \ll 41000
$$

We can see that there is about an order of magnitude difference, and this difference will only get larger as we add convection and thre-dimensionality to the picture. Thus, if is OK to drop $\mathbf{L}\mathbf{D}^{-1}\mathbf{U}$. 

Therefore, return to Eq.(\ref{eq:lusgs-with-additional-term}), which was given as:

$$
\mathbf{A} = (\mathbf{D} + \mathbf{L} + \mathbf{U}) + \mathbf{L}\mathbf{D}^{-1}\mathbf{U}
$$

It is reasonable to write it as:

$$
\mathbf{A} = \mathbf{D} + \mathbf{L} + \mathbf{U}
$$

This might seem that we haven't really gained anything, but remember, our approximate LU decomposition was given by $(\mathbf{D} + \mathbf{L}) \mathbf{D}^{-1} (\mathbf{D} + \mathbf{U})$, and we can say now that this is approximateloy equal to $\mathbf{D} + \mathbf{L} + \mathbf{U}$, i.e.:

$$
\mathbf{D} + \mathbf{L} + \mathbf{U} \approx (\mathbf{D} + \mathbf{L}) \mathbf{D}^{-1} (\mathbf{D} + \mathbf{U})
$$

And this is what is important. We have found an approximation to $\tilde{\mathbf{L}} \, \tilde{\mathbf{U}}$ (which is $(\mathbf{D} + \mathbf{L}) \mathbf{D}^{-1} (\mathbf{D} + \mathbf{U})$), that is formed of only $\mathbf{D}$, $\mathbf{L}$, and $\mathbf{U}$, and approximately equal to $\mathbf{D} + \mathbf{L} + \mathbf{U}$.

OK, we are almost done, so let's write out the final form. What we have done up until this point is to show that:

$$
\mathbf{A} \approx (\mathbf{D} + \mathbf{L}) \mathbf{D}^{-1} (\mathbf{D} + \mathbf{U})
$$

Therefore, if we are solving $\mathbf{Ax}=\mathbf{b}$, we can write this as:

$$
(\mathbf{D} + \mathbf{L}) \mathbf{D}^{-1} (\mathbf{D} + \mathbf{U}) \mathbf{x}=\mathbf{b}
$$

We split this now into two steps. The first one becomes:

$$
(\mathbf{D} + \mathbf{L}) \mathbf{y} = \mathbf{b}
$$

Here, $\mathbf{y}$ is just some intermediate solution array. This is just our forward Gauss-Seidel method. We haven't used $\mathbf{D}^{-1} (\mathbf{D} + \mathbf{U})$ yet, so we write the second step as:

$$
(\mathbf{D} + \mathbf{U})\mathbf{x} = \mathbf{D}\mathbf{y}
$$

This is also the Gauss-Seidel method, but using a backwards substitution (we have the upper triangular matrix $\mathbf{U}$ on the left-hand side, not the lower triangular matrix). Since we do both a forward and backward sweep, and we are using the Gauss-Seidel method here, we have found a solution using the symemtric Gauss-Seidel method.

As I said at the beginning of this section, the LU-SGS method is pointless, we already have the symemtric Gauss-Seidel method, and both are algebraically the same. But, for you as a CFD developer, it means that if you already have the symmetric Gauss-Seidel method implemented, you can now claim that you also have implmeneted the LU-SGS method. You have added a feature without having to code anything.

Unethical? Sure, but most people don't know the difference between the LU-SGS and symmetric Gauss-Seidel method.

## Krylov subspace methods

At this point, we are going to come back to our friend Krylov. We have abandoned him for quite some time and pretended that linear algebra is a field dominated by Gauss (who, last time I checked, is no longer research active). But, it is time now to look at linear system of equation solvers that are somewhat different from what we discussed thus far.

Krylov subspace methods are the de-facto standard for solving linear system of equations in CFD. While you still may find Jacobi, Gauss-Seidel (and its evil twin brother LU-SGS) every now and then (in not yet production ready or academic codes), Krylov subspace methods are the one you want to use when transitioning from simple channel flows or shock tube problems to real-life applications.

I'll spend some time on these method as there are important concepts we need to understand to make sense of them, and to see why these are Krylov subspace methods in the first place. Let's start with the first and, perhaps, most influential method of them all; the Conjugate gradient method.

### The Conjugate Gradient (CG) method 

I find the Conjugate Gradient (CG) method utterly fascinating. So far, we have tried to solve the linear system of equations $\mathbf{Ax}=\mathbf{b}$ directly, by cleverly applying the Gaussian elimination. You will be glad to hear that we are not going to make use of Gauss again. No, the CG method completelty reformulates the problem into an entirely different problem which we then solve.

We need to do a bit of derivation to get to the CG method, and so, we will take our time to develop the CG method. I will break this into 3 separate sections; the reformulation problem, the method of Steepest Descent, and then, the method of conjugate gradients, at last.

#### Reformulating the linear system of equations problem into a minimisation problem

The key idea is to reformulate the the linear system of equations to:

$$
\mathbf{Ax} - \mathbf{b} = 0
$$

We simply write the right-hand side onto the left-side. So far, so good. The idea in the CG method is to replace the left-hand side with something else which is easier to solve. Well, the right-hand side is 0, so we need to think of an equation that we can write on the right-hand side that evaluates to zero.

What about the gradient of a function $f$? If we want to find the minimum or maximum of a function, we take its derivative and set it to zero, right? So, how about we think of a function $f$ and take the derivative of that? Then, if we minimise/maximise the function (that is, we look for its min or max value), we know that the derivative at that point must be zero.

If the derivative of that function is zero, then we can set it equal to $\mathbf{Ax} - \mathbf{b}$. So, let's do that. We write:

$$
\mathbf{Ax} - \mathbf{b} = \nabla f
\tag{eq:minimisation-problem-Ax-b}
$$

Of course, we can't just take any arbitrary function $f$. We need to constrain it based on $\mathbf{A}$, $\mathbf{x}$, and $\mathbf{b}$, i.e. $f$ should be a function of these quantities only.

To do that, we can take Eq.(\ref{eq:minimisation-problem-Ax-b}) and integrate it. This results in:

$$
\int \mathbf{Ax} \mathrm{d}\mathbf{x} - \int \mathbf{b} \mathrm{d}\mathbf{x} = \int\nabla f = f
$$

Indeed, using the integration will give us $f$. So, let's do that. To follow the derivation, it is easiest to do that based on an example. For a simple 2-by-2 problem, we have:

$$
\begin{bmatrix}
\frac{\partial f}{\partial x_1} \\[1em]
\frac{\partial f}{\partial x_2}
\end{bmatrix}=
\begin{bmatrix}
a_{11} & a_{12} \\[1em]
a_{21} & a_{22}
\end{bmatrix}
\begin{bmatrix}
x_{1} \\[1em]
x_{2}
\end{bmatrix} -
\begin{bmatrix}
b_{1} \\[1em]
b_{2}
\end{bmatrix}
\tag{eq:2-by-2-problem-CG}
$$

We will do the derivation for each term independently and then add the results together at the end. So, the first term we will look at is $\int \mathbf{Ax} \mathrm{d}\mathbf{x}$. Let us first write out Eq.(\ref{eq:2-by-2-problem-CG}) for $\mathbf{Ax}$ only (i.e. ignoring $\mathbf{b}$). This gives us:

$$
\nabla f = \mathbf{Ax}\\[1em]
\begin{bmatrix}
\frac{\partial f}{\partial x_1} \\[1em]
\frac{\partial f}{\partial x_2}
\end{bmatrix}=
\begin{bmatrix}
a_{11} & a_{12} \\[1em]
a_{21} & a_{22}
\end{bmatrix}
\begin{bmatrix}
x_{1} \\[1em]
x_{2}
\end{bmatrix}\\[1em]
\begin{bmatrix}
\frac{\partial f}{\partial x_1} \\[1em]
\frac{\partial f}{\partial x_2}
\end{bmatrix}=
\begin{bmatrix}
a_{11}x_1 + a_{12}x_2 \\[1em]
a_{21}x_1 + a_{22}x_2
\end{bmatrix}
\tag{eq:2-by-2-problem-CG-Ax}
$$

From the first row, we can write:

$$
\frac{\partial f}{\partial x_1} = a_{11}x_1 + a_{12}x_2
$$

Now we integrate with respect to $x_1$ and we get:

$$
\int \frac{\partial f}{\partial x_1} \mathrm{d}x_1 = \int a_{11}x_1 \mathrm{d}x_1 + \int a_{12}x_2 \mathrm{d}x_1\\[1em]
f = \frac{a_{11}}{2}x_1^2 + a_{12}x_2 x_1 + g(x_2)
\tag{eq:2-by-2-problem-CG-Ax-row1}
$$

There are two things to take note here. First, $f$ is a scalar function. From Eq.(\ref{eq:2-by-2-problem-CG-Ax}), we can see that we have derivatives of both $\partial f/\partial x_1$ and $\partial f/\partial x_2$, thus, the function $f$ must depend on $f(x_1,x_2)$. Therefore, the integration constant $g(x_2)$ we gain after the integration depends on $x_2$, even though we integrated with respect to $x_1$.

The integration constant $g$ is constant with respect to $x_1$, but, because $f$ depends also on $x_2$, $g$ may not be constant with respect to $x_2$, and so, we write this as $g(x_2)$. If we had a $n$-by-$n$ problem instead, then we would have $g(x_2, x_3, ..., x_n)$.

Since $f$ is a scalar function and depends on both $x_1$ and $x_2$, i.e. its derivative in both direction has to be zero, we can differentiate Eq.(\ref{eq:2-by-2-problem-CG-Ax-row1}) now with respect to $x_2$. tHIS GIVES US:

$$
\frac{\partial f}{\partial x_2} = \frac{\partial}{\partial x_2}\left(\frac{a_{11}}{2}x_1^2 + a_{12}x_2 x_1 + g(x_2)\right)\\[1em]
\frac{\partial f}{\partial x_2} = a_{12}x_1 + g'(x_2)
$$

The first term disappears as it only depends on $x_1$, which is treated as a constant when we differentiate with respect to $x_2$. The derivative of a constant is zero and so it vanishes. The second term depends on both $x_1$ and $x_2$. The derivative of $x_2$ with respect to $x_2$ is $\partial x_2/\partial x_2 = 1$ and so we retain the second term. We also have to differentiate the constant $g(x_2)$, but because we don't know its exact value yet, we simply write this as $\partial g(x_2)/\partial x_2=g'(x_2)$.

Now we go back to Eq.(\ref{eq:2-by-2-problem-CG-Ax}) and evaluate the second row. This gives us:

$$
\frac{\partial f}{\partial x_2} = a_{21}x_1 + a_{22}x_2
$$

Now we have to equations for $\partial f/\partial x_2$ and so we can set them equal and get:

$$
a_{12}x_1 + g'(x_2) = a_{21}x_1 + a_{22}x_2
\tag{eq:2-by-2-problem-CG-Ax-row1-row2-eq}
$$

**Now we make a crucial assumption. We say that our matrix is symmetric.** If that is the case, then $a_{12}=a_{21}$. However, this also means that whatever follows next is only applicable to symmetric matrices (and, symmetric matrices arise if we use numerical schemes with symmetric numerical schemes). For example, the central scheme for second-order derivatives is:

$$
\frac{\partial^2 \phi}{\partial x^2} \approx \frac{\phi_{i+1} - 2\phi_i + \phi_{i-1}}{(\Delta x)^2}
$$

Here, since both $\phi_{i+1}$ and $\phi_{i-1}$ appear with the same factor infront of them ($+1$, which we have dropped here), this numerical scheme will lead to a symemtric matrix when we discretise the equation and form the coefficient matrix $\mathbf{A}$. On the other hand, a non-symmetric example would be the first-order upwind scheme, where we have:

$$
u\frac{\partial \phi}{\partial x}\approx
\begin{cases}
\frac{\phi_i - \phi_{i-1}}{\Delta x} & \text{if } u \ge 0\\[1em]
\frac{\phi_{i+1} - \phi_{i}}{\Delta x} & \text{if } u\lt 0
\end{cases}
$$

Let's say we have $u \ge 0$, then we have $\partial \phi/\partial x \approx \phi_i - \phi_{i-1}/\Delta x$. In this case, we have a factor of $-1$ in front of $\phi_{i-1}$ and a factor of $0$ in front of $\phi_{i+1}$. This means we have a non-symmetric stencil and this would lead to a non-symmetric coefficient matrix.

To summarise, if we restrict ourselves to symmetric matrices only, i.e. we say that $a_{12}=a_{21}$, then we can't use non-symmetric numerical schemes (stencils). As it turns out, this restriction isn't as bad as it first seems, and people have later lifted the restriction and further developed the CG method so that we can also apply this to non-symmetric matrices. But we are getting ahead of ourselves.

Since we said that the coefficient matrix is symmetric, we can write $a_{12}=a_{21}$ as already established. Therefore, from Eq.(\ref{eq:2-by-2-problem-CG-Ax-row1-row2-eq}), we have:

$$
a_{12}x_1 + g'(x_2) = a_{21}x_1 + a_{22}x_2\\[1em]
a_{12}x_1 + g'(x_2) = a_{12}x_1 + a_{22}x_2
$$

Since we now have $a_{12}x_1$ on both sides, we can subtract it from the equation and get:

$$
g'(x_2) = a_{22}x_2
$$

We can integrate this expression now and get:

$$
\int g'(x_2)\mathrm{d}x_2 = \int a_{22}x_2 \mathrm{d}x_2\\[1em]
g(x_2) = \frac{a_{22}}{2}x_2^2
$$

Returning to Eq.(\ref{eq:2-by-2-problem-CG-Ax-row1}), we can insert this integration for $g(x_2)$ and obtain:

$$
f = \frac{a_{11}}{2}x_1^2 + a_{12}x_2 x_1 + \frac{a_{22}}{2}x_2^2
\tag{eq:result-for-f-for-Ax-CG}
$$

Let's look at the second term $a_{12}x_2 x_1$. We want to write this as a matrix vector multiplication in matrix notation, i.e. using $\mathbf{A}$ and $\mathbf{x}$ only. Let's remember that $f$ is a scalar function, so any multiplication of $\mathbf{A}$ and $\mathbf{x}$ needs to result in a scalar. Well, if we simply write $\mathbf{Axx}$, since we have $a_{12}x_2 x_1$, let's see what happens:

$$
\begin{bmatrix}
a_{11} & a_{12} \\[1em]
a_{21} & a_{22}
\end{bmatrix}
\begin{bmatrix}
x_{1} \\[1em]
x_{2}
\end{bmatrix}
\begin{bmatrix}
x_{1} \\[1em]
x_{2}
\end{bmatrix}
$$

Multiplying out the first matrix vector multiplication results in:

$$
\begin{bmatrix}
a_{11}x_1 + a_{12}x_2 \\[1em]
a_{21}x_1 + a_{22}x_2
\end{bmatrix}
\begin{bmatrix}
x_{1} \\[1em]
x_{2}
\end{bmatrix}
$$

Now we have the multiplication of two column vectors, and that is not defined, so we can't replace $a_{12}x_2 x_1$ with $\mathbf{Axx}$. However, we can multiply a column and a row vector together. So, let's see what happens if we rewrite our expression as $\mathbf{Axx^T}$. Skipping the $\mathbf{Ax}$ multiplication step which we have already performed, and evaluation $(\mathbf{Ax})\mathbf{x}^T$ directy, we have:

$$
\begin{bmatrix}
a_{11}x_1 + a_{12}x_2 \\[1em]
a_{21}x_1 + a_{22}x_2
\end{bmatrix}
\begin{bmatrix}
x_{1} & x_{2}
\end{bmatrix}=
\begin{bmatrix}
(a_{11}x_1 + a_{12}x_2)x_1 & (a_{11}x_1 + a_{12}x_2)x_2 \\[1em]
(a_{21}x_1 + a_{22}x_2)x_1 & (a_{21}x_1 + a_{22}x_2)x_2
\end{bmatrix}
$$

Well, while this operation is mathematically defined, what we have ended up doing is to convert our column vector from the product $\mathbf{Ax}$ back into a matrix by multiplying with the transpose (row vector) of $\mathbf{x}^T$. However, we said that $f$ is a scalar function. So, multiplying a column vector with a row vector produces a matrix, but, multiplying a row vector with a column vector produces a scalar!

So, instead of evaluating $\mathbf{(Ax)x^T}$, we need to bring $\mathbf{x}^T$ to the front, i.e. $\mathbf{x}^T(\mathbf{Ax})$. Then, $\mathbf{x}^T$ is a row vector and $(\mathbf{Ax})$ is a column vector. Let's write this out:

$$
\begin{bmatrix}
x_{1} & x_{2}
\end{bmatrix}
\begin{bmatrix}
a_{11}x_1 + a_{12}x_2 \\[1em]
a_{21}x_1 + a_{22}x_2
\end{bmatrix}=(a_{11}x_1 + a_{12}x_2)x_1 + (a_{21}x_1 + a_{22}x_2)x_2
$$

This can be written out explicitly as:

$$
\mathbf{x}^T\mathbf{Ax} = a_{11}x_1^2 + a_{12}x_2x_1 + a_{21}x_1x_2 + a_{22}x_2^2
$$

Since we have $a_{12}=a_{21}$, we can also write this as:

$$
\mathbf{x}^T\mathbf{Ax} = a_{11}x_1^2 + a_{12}x_2x_1 + a_{12}x_1x_2 + a_{22}x_2^2\\[1em]
\mathbf{x}^T\mathbf{Ax} = a_{11}x_1^2 + 2a_{12}x_2x_1 + a_{22}x_2^2 
$$

Let us now multiply both sides by $1/2$. This gives us:

$$
\frac{1}{2}\mathbf{x}^T\mathbf{Ax} = \frac{a_{11}}{2}x_1^2 + a_{12}x_2x_1 + \frac{a_{22}}{2}x_2^2 
$$

Now, compare that against the result we obtained in Eq.(\ref{eq:result-for-f-for-Ax-CG}), which was:

$$
f = \frac{a_{11}}{2}x_1^2 + a_{12}x_2 x_1 + \frac{a_{22}}{2}x_2^2
$$

We can see, these expressions are identical. Therefore, we can write that:

$$
\int \mathbf{Ax}\mathrm{d}\mathbf{x} = \frac{1}{2}\mathbf{x}^T(\mathbf{Ax})
\tag{eq:cg-f-term-1}
$$

Remember what we have done, or tried to achieve here. We first set $\mathbf{Ax} - \mathbf{b} = 0$, and then said that the right-hand side is equal to the gradient of $f$, or $\nabla f$. By integration, we get rid of the gradient operator and then we need to integrate each term individually, and we have just found the first term. Now, let us move on to the second term, this is given as:

$$
\mathbf{b} = \nabla f 
$$

We can also write this as:

$$
\begin{bmatrix}
b_1 \\[1em]
b_2 \\[1em]
\vdots \\[1em]
b_n \\[1em]
\end{bmatrix} =
\begin{bmatrix}
\frac{\partial f}{\partial x_1} \\[1em]
\frac{\partial f}{\partial x_2} \\[1em]
\vdots \\[1em]
\frac{\partial f}{\partial x_n} \\[1em]
\end{bmatrix}
$$

Or, in index notation:

$$
b_i = \frac{\partial f}{\partial x_i}
$$

We can integrate this exprerssion again and get:

$$
\int b_i\mathrm{d}x_i = f
$$

Here, $b_i$ is a constant, and so, we only need to evaluate $\int \mathrm{d}x_i$. If $\mathbf{x}$ consists of only a single entry, then we have $\int\mathrm{d}x_i = \int\mathrm{d}x = x + c$, however, we have as many entries in $\mathbf{x}$ as we have entries in $\mathbf{b}$. Thus, we need to integrate over each entry, leading to:

$$
\int b_i\mathrm{d}x_i = f\\[1em]
(b_1 x_1 + b_2 x_2 + ... + b_n x_n) + c = f\\[1em]
b_i x_i + c = f
$$

Now, we need to transform this index notation again into a matrix notation using $\mathbf{b}$ and $\mathbf{x}$. What we have is the dot product, or scalar product, here, and so we could write:

$$
\mathbf{b}\cdot\mathbf{x} + c= f
$$

However, the literature tends to write it as:

$$
\mathbf{x}^T\mathbf{b} + c= f
\tag{eq:cg-f-term-2}
$$

Both $\mathbf{x}$ and $\mathbf{b}$ are column vectors and, when we evaluated the integral of $\mathbf{Ax}$, we already saw that a row vector multiplied by a column vector is a scalar (hence, we brint the transpose of $\mathbf{x}$ to the left, which makes it a row vector). In fact, for two generic column vectors $\mathbf{a}$ and $\mathbf{b}$, we have the following equality:

$$
\mathbf{a}^T\mathbf{b}=\mathbf{a}\cdot\mathbf{b}
$$

We can show this by writing out the vectors:

$$
\begin{bmatrix}
a_1 & a_2 & ... & a_n
\end{bmatrix}
\begin{bmatrix}
b_1 \\
b_2 \\
\vdots \\
b_n
\end{bmatrix}=
\begin{bmatrix}
a_1 \\
a_2 \\
\vdots \\
a_n
\end{bmatrix}\cdot
\begin{bmatrix}
b_1 \\
b_2 \\
\vdots \\
b_n
\end{bmatrix}
$$

If we now carry out the vector-vector multiplication on the left-hand side, and the dot product on the right-hand side, we have:

$$
a_1 b_1 + a_2 b_2 + ... + a_n b_n = a_1 b_1 + a_2 b_2 + ... + a_n b_n
$$

As we can see, both the left-hand side and right-hand side are identical (by virtue of ```Ctrl+c``` and ```Ctrl+v```).

Returning to Eq.(\ref{eq:cg-f-term-2}), we ses that the integration over $\mathbf{b}$ results in $\mathbf{x}^T\mathbf{b} + c$. We don't know the integration constant, so we just leave it as the letter $c$. We can now combine Eq.(\ref{eq:cg-f-term-1}) and Eq.(\ref{eq:cg-f-term-2}) and insert them into our minimisation problem:

$$
\mathbf{Ax} - \mathbf{b} = 0 \\[1em]
\mathbf{Ax} - \mathbf{b} = \nabla f \\[1em]
\int\mathbf{Ax}\mathrm{d}\mathbf{x} - \int\mathbf{b}\mathrm{d}\mathbf{x} = f \\[1em]
\frac{1}{2}\mathbf{x}^T(\mathbf{Ax}) - \mathbf{x}^T\mathbf{b} + c = f
$$

Or, swapping the left and right side:

$$
f = \frac{1}{2}\mathbf{x}^T(\mathbf{Ax}) - \mathbf{x}^T\mathbf{b} + c
\tag{eq:cg-minimsation-function-f}
$$

While preparing this section, I have gone through my CFD textbooks, as well as other collected material (lecture scripts, lecture presentations, etc.) that deal with the Conjugate Gradient method. Here are three quotes I pulled from the literature:

- *Without derivation, introduce a function $f$:*
- *With a little bit of tedious math, one can apply Equiation 5 to Equation 3, adn derive*
- *Through mathematical manipulations the gradient [of $f$] is*

If they had been more honest, they all could have written:

*I personally don't know how to derive f and I can't be bothered to look it up or even attempt to derive it myself, so, without proof, here is f from a source I found. Let's hope that the authors of the source I am consulting did at least derive it themselves and that they did not just copy it from somewhere else as well.*

Ah yes, as we all know, *hope* is the foundation of scientific progress. I have said it already in the article on deriving the Navier-Stokes equations, but I have issues with textbooks that don't derive their equations. A textbook on CFD that simply lists equations is not a textbook, and yet, most of the books on my shelve resort to this level of lazyness.

I digress.

Let's step back and see what we have done. We have a linear system of equations $\mathbf{Ax}=\mathbf{b}$. Instead of solving this, we solve Eq.(\ref{eq:cg-minimsation-function-f}). Take a look at that equation, it still has our coefficient matrix $\mathbf{A}$, which does not change in this equation. But, we also have our unknown vector $\mathbf{x}$.

So, our goal is to find a vector $\mathbf{x}$ for which $f$ is minimised (i.e. $f$ has the smallest value). Whichever vector $\mathbf{x}$ gives me the lowest function value of $f$ also represents a solution for our linear system of equations, i.e. $\mathbf{Ax}=\mathbf{b}$. So, instead of doing any Gaussian elimination, or even inverting any matrix, we simply solve a minimisation problem.

Now, I have stated a few times now that we minimise $f$, not maximise it. Why is that? In the beginning of this section, I stated that the system of equations can be written as $\mathbf{Ax} - \mathbf{b} = 0 = \nabla f$. Here, $\nabla f$ would evaluate to zero for both a local minimum or maximum (or even a saddle point). So, why do we restrict ourselves to a minimisation?

This is where we make yet one more restriction on our matrix. We already said that the matrix is symmetric, that is $\mathbf{A}^T=\mathbf{A}$, and now we also say that it has to be positive definite. Remember that property? Positive definite matrices have the following property:

$$
\mathbf{x}^T\mathbf{Ax}\gt 0,\qquad \mathbf{x}\ne\mathbf{0}
$$

If this where to be less than zero, we had a negative definite matrix, and then we would have to maximise $f$, instead of minimising it. How can we show that? Well, let's return to Eq.(\ref{eq:cg-minimsation-function-f}). Instead of treating $\mathbf{A}$ like a matrix, let's say it is a 1-by-1 matrix, so it reduces to $\mathbf{A}=a$. Then, our vectors $\mathbf{x}$ and $\mathbf{b}$ simply become scalar values and can be written as $x$ and $b$. We have:

$$
f = \frac{1}{2}\mathbf{x}^T(\mathbf{Ax}) - \mathbf{x}^T\mathbf{b} + c\\[1em]
f = \frac{1}{2}xax - xb + c\\[1em]
f = \frac{1}{2}ax^2 - bx + c
$$

Since we require now $\mathbf{A}$ to be positive definite, we said that we have $\mathbf{x}^T\mathbf{Ax}\gt 0$. For our 1-by-1 matrix, we have $ax^2\gt 0$. Since $x$ is squared, it is always positive. So, $a$ has to be positive. If the leading term in a quadratic function is positive, then we get a parabola that opens up to the top. If it is negative, then it opens to the bottom, as shown in the following example:

<!-- wp:image {"width":"400px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\quadratic_function.webp" alt="This plot shows two quadratic functions, one opening to the top with a being greater than zero, and one opening to the bottom, with a being less than zero." class="wp-image-5550" style="width:400px"/><figcaption class="wp-element-caption">Figure reproduced from <a href="https://www.thethinkacademy.com/blog/how-to-solve-quadratic-equations-forms-and-methods/" target="_blank" rel="noopener" title="">Think academy</a>.</figcaption></figure>
<!-- /wp:image -->

Thus, we always have $a\gt 0$, or $\mathbf{x}^T\mathbf{Ax}\gt 0$ if we restrict ourselves to positive definite matrices, and these typically appear anyways as part of our discretisation procedure. Therefore, we will always have a parabola that opens to the top (or a [hypersurface](https://en.wikipedia.org/wiki/Hypersurface) in higher-dimensional space). Regardless of what type of abstract surface we are looking at, it will always open up to the top, and therefore, it will always have a minimim due to the positive definite restriction.

OK, so we have now derived all the math we need to understand why we can transform a linear system of equations into an equivalent minimisation problem. However, before we crack on with yet more linear algebra, I thought it may be good to look at a problem. I am going to purposefully select a simple problem with a 2-by-2 matrix. This doesn't have any relevance for CFD, unfortunately, but it does mean that we can look at the solution in a 2D space and visually see the equivalence between solving $\mathbf{Ax}=\mathbf{b}$ and minimising $f$.

The example matrix I will be using is:

$$
\mathbf{A}=
\begin{bmatrix}
2 & 1 \\
1 & 2
\end{bmatrix}
$$

The right-hand side vector I'll be using is:

$$
\mathbf{b}=
\begin{bmatrix}
5 \\
4
\end{bmatrix}
$$

Since this example is so simple, we can compute $\mathbf{x}$ directly if we wanted from $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$. This mean we have to form the inverse of $\mathbf{A}^{-1}$, which we have avoided thus far like the pest. For larger systems we don't want to invert $\mathbf{A}$, but for a 2-by-2 system, we should be OK. The inverse of a 2-by-2 matrix can be written as:

$$
\mathbf{A} =
\begin{bmatrix}
a & b \\
c & d
\end{bmatrix}\qquad
\mathbf{A}^{-1} =
\frac{1}{\det({\mathbf{A}})}
\begin{bmatrix}
d & -b \\
-c & a
\end{bmatrix} =
\frac{1}{ad-bc}
\begin{bmatrix}
d & -b \\
-c & a
\end{bmatrix}
$$

Inserting values, we get:

$$
\mathbf{A}^{-1} =
\frac{1}{2\cdot 2-1\cdot 1}
\begin{bmatrix}
2 & -1 \\
-1 & 2
\end{bmatrix} =
\frac{1}{3}
\begin{bmatrix}
2 & -1 \\
-1 & 2
\end{bmatrix}
$$

We can now compute the solution for $\mathbf{x}$ as:

$$
\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}=
\frac{1}{3}
\begin{bmatrix}
2 & -1 \\
-1 & 2
\end{bmatrix}
\begin{bmatrix}
5 \\
4
\end{bmatrix}=
\frac{1}{3}
\begin{bmatrix}
2\cdot 5 + (-1\cdot 4) \\
(-1\cdot 5) + 2\cdot 4
\end{bmatrix}=
\frac{1}{3}
\begin{bmatrix}
10 + (-4) \\
(-5) + 8
\end{bmatrix}=
\frac{1}{3}
\begin{bmatrix}
6 \\
3
\end{bmatrix}=
\begin{bmatrix}
2 \\
1
\end{bmatrix}
$$

Thus, we now know that the solution should be $\mathbf{x}=[2,1]^T$. With that knowledge, let's evaluate the function of $f$ and see what happens. Remember, $f$ was given as:

$$
f = \frac{1}{2}\mathbf{x}^T(\mathbf{Ax}) - \mathbf{x}^T\mathbf{b} + c
$$

Since we want to minimise this function, the constant $c$ is just going to raise or lower the minimum value, but it will not affect at which location $\mathbf{x}$ this occurs. So, we typically set this value to $c=0$ and then we don't have to worry about it anymore (in-fact, we can remove it completely; it won't affect the solution we find).

The code below evaluates this function in the interval of $-5\le x \le 5$ and $-5 \le y \le 5$. This interval contains our solution $\mathbf{x}=[2,1]^T$, i.e. the solution is at $x=2$ and $y=1$ and so, we should find a minimum of $f$ at that point. Here is the (python) code with comments:

<!-- wp:kevinbatdorf/code-block-pro {"code":"import matplotlib.pyplot as plt\nimport numpy as np\n\n# create 2-by-2 matrix\nA = np.array([[2, 1], [1, 2]])\n\n# right-hand side vector\nb = np.array([5, 4])\n\n# solving for x using x = A^{-1}b, which is cheap for a 2-by-2 matrix\nxsol = np.linalg.solve(A, b) \n\n# solution at [2 1]\nprint('solution at ', xsol)\n\n# let's create a search space for f between -5 and 5 for both x and y\n# we know that the solution is at [2 1], which is within this interval\nxc = np.linspace(-5, 5, 100)\nyc = np.linspace(-5, 5, 100)\n\n# creating an empty array for f, where we store the solutions for f so we can plot it later\nf = np.zeros((len(xc), len(yc)))\n\n# the constant is not impoortant, it will just shift our plane up or down\n# it won't have an effect so we can set it to 0\nc = 0\n\n# compute for all values of x and y\nfor i in range(0, len(xc)):\n    for j in range(0, len(yc)):\n        xvec = np.array([xc[i], yc[j]])\n\n        # change i,j to j,i indices here due to difference in how indices\n        # are interpreted in numpy and matplot lib\n        f[j][i] = 0.5 * xvec.T @ A @ xvec - b.T @ xvec + c\n\n# create grid for surface\nX, Y = np.meshgrid(xc, yc)\n\n# create figure with 1x2 layout\nfig = plt.figure(figsize=(10,5))\n\n# \u002d\u002d\u002d\u002d LEFT: contour plot \u002d\u002d\u002d\u002d\nax1 = fig.add_subplot(1, 2, 1)\n\ncont = ax1.contourf(xc, yc, f, 50, cmap='viridis')\n\n# plot solution\nax1.plot(xsol[0], xsol[1], 'rx')\n\n# dashed helper lines\nax1.plot([2, 2], [-5, 0.7], 'r\u002d\u002d')\nax1.plot([-5, 1.7], [1, 1], 'r\u002d\u002d')\n\nax1.set_xlabel('x')\nax1.set_ylabel('y')\n\nfig.colorbar(cont, ax=ax1)\n\n# \u002d\u002d\u002d\u002d RIGHT: 3D surface \u002d\u002d\u002d\u002d\nax2 = fig.add_subplot(1, 2, 2, projection='3d')\n\nax2.plot_surface(X, Y, f, cmap='viridis')\n\nax2.set_xlabel('x')\nax2.set_ylabel('y')\nax2.set_zlabel('f')\n\nax2.view_init(elev=30, azim=-30)\n\nplt.tight_layout()\nplt.show()","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003eimport\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e matplotlib.pyplot \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003eas\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e plt\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003eimport\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e numpy \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003eas\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e np\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# create 2-by-2 matrix\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eA = np.array([\u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;, \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;])\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# right-hand side vector\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eb = np.array(\u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e4\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# solving for x using x = A^{-1}b, which is cheap for a 2-by-2 matrix\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003exsol = np.linalg.solve(A, b) \u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# solution at \u0026#91;2 1\u0026#93;\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003eprint\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;solution at \u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, xsol)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# let\u0026#39;s create a search space for f between -5 and 5 for both x and y\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# we know that the solution is at \u0026#91;2 1\u0026#93;, which is within this interval\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003exc = np.linspace(-\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e100\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eyc = np.linspace(-\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e100\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# creating an empty array for f, where we store the solutions for f so we can plot it later\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003ef = np.zeros((\u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003elen\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(xc), \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003elen\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(yc)))\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# the constant is not impoortant, it will just shift our plane up or down\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# it won\u0026#39;t have an effect so we can set it to 0\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003ec = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# compute for all values of x and y\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e i \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003ein\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003erange\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003elen\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(xc)):\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e    \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003efor\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e j \u003c/span\u003e\u003cspan style=\u0022color: #C586C0\u0022\u003ein\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003erange\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #DCDCAA\u0022\u003elen\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e(yc)):\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        xvec = np.array([xc\u0026#91;i\u0026#93;, yc\u0026#91;j\u0026#93;])\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# change i,j to j,i indices here due to difference in how indices\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        \u003c/span\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# are interpreted in numpy and matplot lib\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e        f\u0026#91;j\u0026#93;\u0026#91;i\u0026#93; = \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0.5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e * xvec.T @ A @ xvec - b.T @ xvec + c\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# create grid for surface\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eX, Y = np.meshgrid(xc, yc)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# create figure with 1x2 layout\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003efig = plt.figure(\u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003efigsize\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e=(\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e10\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e,\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e))\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# \u002d\u002d\u002d\u002d LEFT: contour plot \u002d\u002d\u002d\u002d\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax1 = fig.add_subplot(\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003econt = ax1.contourf(xc, yc, f, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e50\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003ecmap\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e=\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;viridis\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# plot solution\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax1.plot(xsol\u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;, xsol\u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;, \u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;rx\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# dashed helper lines\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax1.plot(\u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;, \u0026#91;-\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e0.7\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;, \u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;r\u002d\u002d\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax1.plot(\u0026#91;-\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e5\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1.7\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;, \u0026#91;\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u0026#93;, \u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;r\u002d\u002d\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax1.set_xlabel(\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;x\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax1.set_ylabel(\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;y\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003efig.colorbar(cont, \u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003eax\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e=ax1)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #6A9955\u0022\u003e# \u002d\u002d\u002d\u002d RIGHT: 3D surface \u002d\u002d\u002d\u002d\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax2 = fig.add_subplot(\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e1\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e2\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003eprojection\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e=\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;3d\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax2.plot_surface(X, Y, f, \u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003ecmap\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e=\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;viridis\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax2.set_xlabel(\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;x\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax2.set_ylabel(\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;y\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax2.set_zlabel(\u003c/span\u003e\u003cspan style=\u0022color: #CE9178\u0022\u003e\u0026#39;f\u0026#39;\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eax2.view_init(\u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003eelev\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e=\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e30\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e, \u003c/span\u003e\u003cspan style=\u0022color: #9CDCFE\u0022\u003eazim\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e=-\u003c/span\u003e\u003cspan style=\u0022color: #B5CEA8\u0022\u003e30\u003c/span\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e)\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eplt.tight_layout()\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eplt.show()\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"python","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":72,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(2 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly>import matplotlib.pyplot as plt
import numpy as np

# create 2-by-2 matrix
A = np.array([&#91;2, 1&#93;, &#91;1, 2&#93;])

# right-hand side vector
b = np.array(&#91;5, 4&#93;)

# solving for x using x = A^{-1}b, which is cheap for a 2-by-2 matrix
xsol = np.linalg.solve(A, b) 

# solution at &#91;2 1&#93;
print('solution at ', xsol)

# let's create a search space for f between -5 and 5 for both x and y
# we know that the solution is at &#91;2 1&#93;, which is within this interval
xc = np.linspace(-5, 5, 100)
yc = np.linspace(-5, 5, 100)

# creating an empty array for f, where we store the solutions for f so we can plot it later
f = np.zeros((len(xc), len(yc)))

# the constant is not impoortant, it will just shift our plane up or down
# it won't have an effect so we can set it to 0
c = 0

# compute for all values of x and y
for i in range(0, len(xc)):
    for j in range(0, len(yc)):
        xvec = np.array([xc&#91;i&#93;, yc&#91;j&#93;])

        # change i,j to j,i indices here due to difference in how indices
        # are interpreted in numpy and matplot lib
        f&#91;j&#93;&#91;i&#93; = 0.5 * xvec.T @ A @ xvec - b.T @ xvec + c

# create grid for surface
X, Y = np.meshgrid(xc, yc)

# create figure with 1x2 layout
fig = plt.figure(figsize=(10,5))

# ---- LEFT: contour plot ----
ax1 = fig.add_subplot(1, 2, 1)

cont = ax1.contourf(xc, yc, f, 50, cmap='viridis')

# plot solution
ax1.plot(xsol&#91;0&#93;, xsol&#91;1&#93;, 'rx')

# dashed helper lines
ax1.plot(&#91;2, 2&#93;, &#91;-5, 0.7&#93;, 'r--')
ax1.plot(&#91;-5, 1.7&#93;, &#91;1, 1&#93;, 'r--')

ax1.set_xlabel('x')
ax1.set_ylabel('y')

fig.colorbar(cont, ax=ax1)

# ---- RIGHT: 3D surface ----
ax2 = fig.add_subplot(1, 2, 2, projection='3d')

ax2.plot_surface(X, Y, f, cmap='viridis')

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('f')

ax2.view_init(elev=30, azim=-30)

plt.tight_layout()
plt.show()</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #C586C0">import</span><span style="color: #D4D4D4"> matplotlib.pyplot </span><span style="color: #C586C0">as</span><span style="color: #D4D4D4"> plt</span></span>
<span class="line"><span style="color: #C586C0">import</span><span style="color: #D4D4D4"> numpy </span><span style="color: #C586C0">as</span><span style="color: #D4D4D4"> np</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># create 2-by-2 matrix</span></span>
<span class="line"><span style="color: #D4D4D4">A = np.array([&#91;</span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">&#93;, &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">&#93;])</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># right-hand side vector</span></span>
<span class="line"><span style="color: #D4D4D4">b = np.array(&#91;</span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">4</span><span style="color: #D4D4D4">&#93;)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># solving for x using x = A^{-1}b, which is cheap for a 2-by-2 matrix</span></span>
<span class="line"><span style="color: #D4D4D4">xsol = np.linalg.solve(A, b) </span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># solution at &#91;2 1&#93;</span></span>
<span class="line"><span style="color: #DCDCAA">print</span><span style="color: #D4D4D4">(</span><span style="color: #CE9178">&#39;solution at &#39;</span><span style="color: #D4D4D4">, xsol)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># let&#39;s create a search space for f between -5 and 5 for both x and y</span></span>
<span class="line"><span style="color: #6A9955"># we know that the solution is at &#91;2 1&#93;, which is within this interval</span></span>
<span class="line"><span style="color: #D4D4D4">xc = np.linspace(-</span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">100</span><span style="color: #D4D4D4">)</span></span>
<span class="line"><span style="color: #D4D4D4">yc = np.linspace(-</span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">100</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># creating an empty array for f, where we store the solutions for f so we can plot it later</span></span>
<span class="line"><span style="color: #D4D4D4">f = np.zeros((</span><span style="color: #DCDCAA">len</span><span style="color: #D4D4D4">(xc), </span><span style="color: #DCDCAA">len</span><span style="color: #D4D4D4">(yc)))</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># the constant is not impoortant, it will just shift our plane up or down</span></span>
<span class="line"><span style="color: #6A9955"># it won&#39;t have an effect so we can set it to 0</span></span>
<span class="line"><span style="color: #D4D4D4">c = </span><span style="color: #B5CEA8">0</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># compute for all values of x and y</span></span>
<span class="line"><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> i </span><span style="color: #C586C0">in</span><span style="color: #D4D4D4"> </span><span style="color: #DCDCAA">range</span><span style="color: #D4D4D4">(</span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #DCDCAA">len</span><span style="color: #D4D4D4">(xc)):</span></span>
<span class="line"><span style="color: #D4D4D4">    </span><span style="color: #C586C0">for</span><span style="color: #D4D4D4"> j </span><span style="color: #C586C0">in</span><span style="color: #D4D4D4"> </span><span style="color: #DCDCAA">range</span><span style="color: #D4D4D4">(</span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">, </span><span style="color: #DCDCAA">len</span><span style="color: #D4D4D4">(yc)):</span></span>
<span class="line"><span style="color: #D4D4D4">        xvec = np.array([xc&#91;i&#93;, yc&#91;j&#93;])</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #6A9955"># change i,j to j,i indices here due to difference in how indices</span></span>
<span class="line"><span style="color: #D4D4D4">        </span><span style="color: #6A9955"># are interpreted in numpy and matplot lib</span></span>
<span class="line"><span style="color: #D4D4D4">        f&#91;j&#93;&#91;i&#93; = </span><span style="color: #B5CEA8">0.5</span><span style="color: #D4D4D4"> * xvec.T @ A @ xvec - b.T @ xvec + c</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># create grid for surface</span></span>
<span class="line"><span style="color: #D4D4D4">X, Y = np.meshgrid(xc, yc)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># create figure with 1x2 layout</span></span>
<span class="line"><span style="color: #D4D4D4">fig = plt.figure(</span><span style="color: #9CDCFE">figsize</span><span style="color: #D4D4D4">=(</span><span style="color: #B5CEA8">10</span><span style="color: #D4D4D4">,</span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">))</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># ---- LEFT: contour plot ----</span></span>
<span class="line"><span style="color: #D4D4D4">ax1 = fig.add_subplot(</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">cont = ax1.contourf(xc, yc, f, </span><span style="color: #B5CEA8">50</span><span style="color: #D4D4D4">, </span><span style="color: #9CDCFE">cmap</span><span style="color: #D4D4D4">=</span><span style="color: #CE9178">&#39;viridis&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># plot solution</span></span>
<span class="line"><span style="color: #D4D4D4">ax1.plot(xsol&#91;</span><span style="color: #B5CEA8">0</span><span style="color: #D4D4D4">&#93;, xsol&#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">&#93;, </span><span style="color: #CE9178">&#39;rx&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># dashed helper lines</span></span>
<span class="line"><span style="color: #D4D4D4">ax1.plot(&#91;</span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">&#93;, &#91;-</span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">0.7</span><span style="color: #D4D4D4">&#93;, </span><span style="color: #CE9178">&#39;r--&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"><span style="color: #D4D4D4">ax1.plot(&#91;-</span><span style="color: #B5CEA8">5</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">1.7</span><span style="color: #D4D4D4">&#93;, &#91;</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">&#93;, </span><span style="color: #CE9178">&#39;r--&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">ax1.set_xlabel(</span><span style="color: #CE9178">&#39;x&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"><span style="color: #D4D4D4">ax1.set_ylabel(</span><span style="color: #CE9178">&#39;y&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">fig.colorbar(cont, </span><span style="color: #9CDCFE">ax</span><span style="color: #D4D4D4">=ax1)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #6A9955"># ---- RIGHT: 3D surface ----</span></span>
<span class="line"><span style="color: #D4D4D4">ax2 = fig.add_subplot(</span><span style="color: #B5CEA8">1</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #B5CEA8">2</span><span style="color: #D4D4D4">, </span><span style="color: #9CDCFE">projection</span><span style="color: #D4D4D4">=</span><span style="color: #CE9178">&#39;3d&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">ax2.plot_surface(X, Y, f, </span><span style="color: #9CDCFE">cmap</span><span style="color: #D4D4D4">=</span><span style="color: #CE9178">&#39;viridis&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">ax2.set_xlabel(</span><span style="color: #CE9178">&#39;x&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"><span style="color: #D4D4D4">ax2.set_ylabel(</span><span style="color: #CE9178">&#39;y&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"><span style="color: #D4D4D4">ax2.set_zlabel(</span><span style="color: #CE9178">&#39;f&#39;</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">ax2.view_init(</span><span style="color: #9CDCFE">elev</span><span style="color: #D4D4D4">=</span><span style="color: #B5CEA8">30</span><span style="color: #D4D4D4">, </span><span style="color: #9CDCFE">azim</span><span style="color: #D4D4D4">=-</span><span style="color: #B5CEA8">30</span><span style="color: #D4D4D4">)</span></span>
<span class="line"></span>
<span class="line"><span style="color: #D4D4D4">plt.tight_layout()</span></span>
<span class="line"><span style="color: #D4D4D4">plt.show()</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

This code plots two things at the end; a contour plot (i.e. a 2D projection of a 3D surface) and the 3D surface itself.

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\f_minimisation_example.png" alt="This plot shows the surface plot of f, where the minimum corresponds to exactly the same location where the solution of x is found." class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

On the left, I show the contours of $f$ between $-5\le x \le 5$ and $-5 \le y \le 5$. I have shown, with a red cross, where the point $x=2$ and $y=1$ is located. These are the coordinates from the solution vector we have compued analytically (the code does also compute this solution numerically as ```xsol = np.linalg.solve(A, b)```). We can see, indeed, this is where the function $f$ appears to have its lowest function values, based on the colorbar.

I have also shown the surface in a 3D plot to the right. It is more difficult to see here that the minimum of this surface is indeed at $x=2$ and $y=1$, but we can at least see that there is a minimum at this location. 

What we can also see is that, indeed, we have a minimum, not a maximum, and so our matrix must be positive definite. We can check that for our solution vector and evaluate

$$
\mathbf{x}^T\mathbf{Ax} =
\begin{bmatrix}
2 & 1
\end{bmatrix}
\begin{bmatrix}
2 & 1\\
1 & 2
\end{bmatrix}
\begin{bmatrix}
2\\
1
\end{bmatrix}=
\begin{bmatrix}
2 & 1
\end{bmatrix}
\begin{bmatrix}
2\cdot 2 + 1\cdot 1\\
1\cdot 2 + 2\cdot 1
\end{bmatrix}=
\begin{bmatrix}
2 & 1
\end{bmatrix}
\begin{bmatrix}
5\\
4
\end{bmatrix}=
2\cdot 5 + 4\cdot 1 = 10 + 4 = 14
$$

Indeed, for this particular point, we have shown that $\mathbf{x}^T\mathbf{Ax}\gt 0$, and you can insert any other vector $\mathbf{x}$, it will always give you a positive value. So, the matrix is not just symmetric ($\mathbf{A}^T=\mathbf{A}$), but it is also positive definite.

Right, so we have spend quite a bit now on showing that we can minimise $f$, instead of solving $\mathbf{Ax}=\mathbf{b}$, so, now we need to talk about algorithms to solve $f$. And this is where we start to get towards the conjugate gradient method. But, before we do that, we will look at an intermediate method, which is known as the method of Steepest Descent.

#### The method of Steepest Descent

The method of Steepest Descents is conceptually speaking very simple. It is just difficult to visualise for anything that doesn't fit into our 3D space. Just looking at the equations doesn't really build up an intuition, and so, we need to approach this from a 3D space initially. The good news is that this is not a limitation, the method generalises to higher dimensions by simple increasing the size of vectors and matrices.

That is, if we have a vector $\mathbf{x}$ with 3 entries, then we can visualise this vector in 3D space. But, if we think of $\mathbf{x}$ as being the solution vector of $\mathbf{Ax}=\mathbf{b}$ arising from a CFD problem, then this means that we have only three entries in our solution vector, e.g. we only have 3 cells in our mesh, and we store the solution inside each cell.

Typically, we have many more entries in $\mathbf{x}$, but if we have a mesh with 1 million cells, so that $\mathbf{x}$ has 1 million entries, we can no longer visualise that, so we have to rely on the maths alone. But, what we derive works just as well in 3D space ($\mathbf{x}$ has 3 entries) as it does for a vector with 1 million entries. The maths is exactly the same, the only difference is the size of the vectors an matrices.

OK, now that we have established that, I want to go back to our function $f$ and now try to minimise that. Let's develop our intuition and see how we can approach this problem.

Initially, I have no idea what the solution for $\mathbf{x}$ is. The whole point is that I minimise $f$, which then tells me what $\mathbf{x}$ is. Thus, to start the algorithm, I have to pick a random value for $\mathbf{x}$. Any value will do. We may pick a zero vector, i.e. one where all entries are zero, or we may set $\mathbf{x}$ equal to a solution computed from a previous timestep/iteration if available.

Once I have picked a starting location, I then need to find a way to update my location $\mathbf{x}$. A simple approach may be the following:

$$
\mathbf{x}^{k+1}=\mathbf{x}^k + \alpha^k\mathbf{r}^k
\tag{eq:line-search}
$$

Here, $\mathbf{x}^{k+1}$ is the updated location after one step. $\mathbf{x}^k$ is the location at the previous step (or my initial guess if we are just starting the algorithm). $\mathbf{r}^k$ is the direction in which we want to go, and $\alpha^k$ tells us by how much we should be going in that direction, i.e. it is just a scalar multiplier.

Now we need to come up with some rules, i.e. how do we obtain the direction vector $\mathbf{r}$, and by how much do we multiply it. I'll propose a stupid algorithm (not to be confused with the [STUPID turbulence model](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/all-you-need-to-know-about-rans-turbulence-modelling-in-one-article/#aioseo-lets-create-our-own-turbulence-model-the-stupid-model)), which I will call "Tom is always right".

No, this is not a political statement, it tells us that we always have to turn right. So, after we have picked an initial direction (any direction, it can be a randomly generated vector), any subsequent vector $\mathbf{r}$ we define has to point to the right. By how much? Well, if we now also generate a random scalar multiplier $\alpha$, we can evaluate Eq.(\ref{eq:line-search}) and get a new value for $\mathbf{x}^{k+1}$.

With this new value for $\mathbf{x}^{k+1}$, we can evaluate $f^{k+1}$ again, and, if this is smaller than the previous value of $f$, i.e. $f^k$ obtained from $\mathbf{x}^k$, then we accept the randomly generated value for $\alpha$ and right-pointing direction vector $\mathbf{r}$. The following figure shows one example of this direction, eventually finding the minimum.

<!-- wp:image {"width":"500px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\line_search_algorithm.png" alt="This figure shows the same function f as seen previously with contours indicating the levels of f in a range of -5 to 5 for both x and y, with x=2 and y=1 being the location of the minimum of f. A random walk shows how one may get to the minimim of f, indicated by points joined together through red dashed lines." class="wp-image-5550" style="width:500px"/></figure>
<!-- /wp:image -->

Now, of course, this is a stupid algorithm (I said as much), and there are many problems with it. First of all, our convergence rate is random. Run this experiment a few times and you will get wildely different results. Secondly, we have to potentially evaluate $f$ many times before we accept a new direction, making the cost per iteration unpredictable. Thirdly, when do we stop? This algorithm will struggle to turn right close to the minimum of $f$, and judging when this minimum is reached is not a trivial task.

The point is, with a simple line search algorithm, as given by Eq.(\ref{eq:line-search}), we can *walk* towards the minimum. In fact, this is exactly what the Conjugate Gradient method will do later (as well as the steepoest decent algorithm). What differentiates both methods is how the determine the direction $\mathbf{r}^k$ in which we go, as well as by how much we go in that direction, that is, the value for the scalar multiplier $\alpha$.

Clearly, we can do better than Tom's stupid algorithm, and thankfully, people with greater intellect before us have solved that issue for us already. This is where the method of Steepest Descent comes in. Let's have a look at a simple parabola $f(x)$, as shown in the following figure:

<!-- wp:image {"width":"500px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\steepest_decent_gradient_direction.png" alt="This figure shows a parabola, which is denoted by f(x), and its derivative, denote by f'(x), at a point to the right of its minimim. The parabola opens to the top. The derivative is positive and shows that going in the negative derivative direction will get us to the minimum." class="wp-image-5550" style="width:500px"/></figure>
<!-- /wp:image -->

We have a random initial location given for $x$ (as indicated by the point at $x=5$), and we want to move to the minimim of $f(x)$. If we look at the derivative of $f(x)$, that is, $f'(x)$, we see that the derivative is positive, i.e. $f'(x) \gt 0$. In fact, with the values provides, we can compute $f'(x)$ as $f'(x)=\Delta y/\Delta x$. For example, evaluating the derivative between $x=4$ and $x=5$, we get:

$$
f'(x)=\frac{\Delta y}{\Delta x}=\frac{1.25-0.25}{5-4} = \frac{1}{1} = 1
$$

So, the derivative is indeed positive. But, if we walked now along the direction of the gradient, we would be going to the right (positive x-direction). This would take us away from the minimum. Instead, we go in the negative direction of the derivative, i.e. $-f'(x)=-1$. By how much are we going in that direction (i.e. what is our scalar multiplier $\alpha$)? Let's simply take the absolute value of the derivative, i.e. $|f'(x)|$. As we get closer to the minimum of $f(x)$, the derivative gets smaller and so our updated will get smaller as well.

We can now use Eq.(\ref{eq:line-search}) again, which becomes:

$$
\mathbf{x}^{k+1}=\mathbf{x}^k + \alpha^k\mathbf{r}^k\\
x^{k+1}=x^k + \alpha^k r^k\\
x^{k+1}=x^k + |f'(x)|\cdot(-f'(x))\\
x^{k+1}=5 + 1\cdot(-1)\\
x^{k+1}=5 -1\\
x^{k+1}=4
$$

So, our new location is $x^{k+1}=4$. We can now evaluate $f(x=4)$ and get a new value for $y$. Then, we can evaluate the derivative again, update $x^{k+1}$ again (where $x^k$ will now be $x^k = 4$, i.e. the value from the previous iteration), and we do that until we have $x^{k+1}\approx x^k$ (i.e. they differ by some small difference only). Eventually, we will get a value of $x^{k+1}\approx 3.5$.

For this simple function depending on only $x$, we could stop here, but we need to now look at the 3D case (or, rather, the case where $f$ depends now on multiple directions). The Steepest Descent algorithm is schematically shown below:

<!-- wp:image {"width":"500px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\steepest_decent.png" alt="This figure shows a surface in 3D space, which has a minimum and opens up like a parabola on all sides to the top. An intersection of the surface with a plane is shown (i.e. a line on the surface), which has a minimum. Another intersection is shown at the minimum, producing another line orthogonal to the first intersection." class="wp-image-5550" style="width:500px"/></figure>
<!-- /wp:image -->

Here, we have again some starting point which may be randomly guessed or given from a previous timestep/iteration, which is indicated by $x^k$. We want to do a line search again, that is, along an arbitrary direction $\mathbf{r}$, we want to go along that direction until we hit a local minimim on the surface (along the direction of $\mathbf{r}$). This will not be the global minimim of $f$, but it will be the minimim along $\mathbf{r}$ starting from $\mathbf{x}^k$, as indicated by the dashed line through $\mathbf{x}^k$.

Now comes a crucial step in the Steepest Descent algorithm: to obtain the next line search direction, we look for the gradient in a direction that is normal (at a 90 degree angle) to the previous search direction. This is shown by the dashed line going from $\mathbf{x}^{k+1}$ to $\mathbf{x}^{k+2}$. We repeat our line search now by finding the local minimum on this new dashed line.

This algorithm converges very quickly to the minimum of $f$, which makes it a very robust algorithm. The question is, how do we find our search direction $\mathbf{r}$ and the scalar multiplier $\alpha$? Well, let's have a look.

To start the algorithm, we need to find a starting point $\mathbf{x}^k$ and a search direction $\mathbf{r}$. The starting point is simple; remember that $\mathbf{x}$ represents our solution vector (e.g. our velocity, pressure, temperature, density, turbulent kinetic energy, etc.), i.e. it is the quantity we are solving for. When we start our simulation, we have to initialise our solution (for example, setting all values the ambient values or values set at an inlet).So, when we start the simulation, $\mathbf{x}$ is just our initial solution.

As the simulation progresses, we set $\mathbf{x}$ equal to the solution we had at the previous iteration, though we could also set $\mathbf{x}$ to anything we want, really. But, setting it to the solution from the previous iteration means that we are probably already quite close to the minimum, and so, we will need fewer iterations to find the minimum.

What about the search direction $\mathbf{r}$? Well, we can construct it from our linear system $\mathbf{Ax}=\mathbf{b}$. We can bring $\mathbf{b}$ onto our left-hand side, and have:

$$
\mathbf{Ax}-\mathbf{b}=0
$$

We saw this before, when we derived the function for $f$. When we did that, we said the right-hand side is equal to the derivative of $f$, i.e. if the derivative is zero, then we have found a local minimum. Thus, we have:

$$
\mathbf{Ax}-\mathbf{b}=\nabla f
$$

When we looked at the method of Steepest Descent, we said that we always want to go in the negative direction of the gradient (as we did in the figure above with the parabola). Therefore, we multiply this equation by $-1$ and have:

$$
-\mathbf{Ax}+\mathbf{b}=-\nabla f
$$

And this is our search direction. We can write this a bit cleaner as:

$$
-\nabla f = -\mathbf{Ax}+\mathbf{b}\\
-\nabla f = \mathbf{b} -\mathbf{Ax}\\
\mathbf{r} = \mathbf{b} -\mathbf{Ax}
\tag{eq:steepest-decent-r}
$$

Ok, so we now know what our starting point $\mathbf{x}^k$ is, we know what the initial search direction is, based on the negative gradient of $f$, but we don't know yet by how much we should go along our search direction $\mathbf{r}$, i.e. we need to find the scalar multiplier $\alpha$. To find that, we need to construct an equation for the search direction update, find a way to get $\alpha$ into that equation, and then solve for it.

In the method of the Steepest Descent, we said that the next search direction should be at an normal angle (90 degrees) to the previous search direction. How can we express that? Well, let's take a simple example in 2D. If my search direction is $\mathbf{r}^k = [1,0]^T$, that is, along the x-direction, which vector would be normal to that? Well, if $\mathbf{r}^k$ is along the x-direction of my coordinate system, a vector normal to that is the y-direction. So, we would expect our next search direction $\mathbf{r}^{k+1}$ to be along the y direction.

It can either point along the positive or negative direction of y, i.e. $\mathbf{r}^{k+1}=[0,1]^T$ and $\mathbf{r}^{k+1}=[0,-1]^T$ are equally possible, it depends on the gradient of $f$ where we evaluate $\mathbf{r}^{k+1}$. When we have two vectors which are normal to each other, as we have in this case, we know that their dot (or scalar) product is zero, that is:

$$
\mathbf{r}^k\cdot\mathbf{r}^{k+1} = 0\\[1em]
\begin{bmatrix}
1\\
0
\end{bmatrix}\cdot
\begin{bmatrix}
0\\
\pm 1
\end{bmatrix}=
1\cdot 0 + 0\cdot(\pm 1) = 0 + 0 = 0
$$

We know that the dot product gives us the angle between two vectors using the cosine function, and, the cosine of 90 degrees is 0, i.e. $\cos (90^\circ)=0$. So, if the dot product evaluates to 0, then the vectors are at a 90 degree angle to each other. We can use this to write our search directions as:

$$
\mathbf{r}^{k+1} \cdot \mathbf{r}^{k}=0
$$

This statement enforces that a new search direction at $k+1$ is going to be normal to our previous search direction at $k$. We have also seen that we can write this dot product in a slightly different form:

$$
\left(\mathbf{r}^{k+1}\right)^T\mathbf{r}^{k}=0
\tag{eq:steepest-decent-search-direction-update}
$$

Here, we writing the search direction at $k+1$ as a row vector, and the search direction at $k$ as a column vector. As we have seen previously, this product is the same as forming the dot product. I'll use this notation from now on, as this is how you will find it in the literature, but know that you can always replace it with the dot product if you wanted.

Remember what we want to do here, we want to find the scalar multiplier $\alpha$ for the current search direction, i.e. we are still at iteration $k$, so we don't think about updating our search direction yet. Therefore, we can use Eq.(\ref{eq:steepest-decent-search-direction-update}) now and see if we can construct an equation for $\alpha$. 

We start by eliminating quantities we don't know. We don't know the search direction $\mathbf{r}^{k+1}$, and so we need to get rid of it. We already saw that the search direction can be computed from Eq.(\ref{eq:steepest-decent-r}). In this equation, we didn't introduce yet indices for $k$ and $k+1$, so let's do this first. We had:

$$
\mathbf{r} = \mathbf{b} -\mathbf{Ax}
$$

If we want to express that as a search direction at $k+1$, then we have

$$
\mathbf{r}^{k+1} = \mathbf{b} -\mathbf{Ax}^{k+1}
$$

If we needed the search direction at $k$, we could have equally set $\mathbf{r}$ and $\mathbf{x}$ at $k$. However, with this equation given now at iteration $k+1$, we can insert that into Eq.(\ref{eq:steepest-decent-search-direction-update}) and get:

$$
\left(\mathbf{b} -\mathbf{Ax}^{k+1}\right)^T\mathbf{r}^{k}=0
\tag{eq:steepest-decent-almost-done}
$$

Hmmm, ok, so we got rid of $\mathbf{r}^{k+1}$, but now we have the updated location at $\mathbf{x}^{k+1}$ in the equation. We don't know yet where it is. The whole point of finding an equation for $\alpha$ is that we can compute what $\mathbf{x}^{k+1}$ will be. So, sorry, $\mathbf{x}^{k+1}$, you have to go as well. If we look back at our line search equation, i.e. Eq.(\ref{eq:line-search}), we had:

$$
\mathbf{x}^{k+1}=\mathbf{x}^k + \alpha^k\mathbf{r}^k
$$

Ah, perfect, we can get rid of $\mathbf{x}^{k+1}$ and, at the same time, replace it not just with an equation for $\alpha$, but all other quantities are known in this expression. Cool, let's insert that then into Eq.(\ref{eq:steepest-decent-almost-done}). This gives us:

$$
\left(\mathbf{b} -\mathbf{A}\left[\mathbf{x}^k + \alpha^k\mathbf{r}^k\right]\right)^T\mathbf{r}^{k}=0
$$

So, we now have an equation where all quantities are known, except for $\alpha$. This means we can solve the equation for $\alpha$ and compute our scalar multiplier. We start by multiplying out the parentheses:

$$
\left(\mathbf{b} -\mathbf{A}\left[\mathbf{x}^k + \alpha^k\mathbf{r}^k\right]\right)^T\mathbf{r}^{k}=0\\[1em]
\left(\mathbf{b} -\mathbf{A}\mathbf{x}^k - \alpha^k\mathbf{A}\mathbf{r}^k\right)^T\mathbf{r}^{k}=0\\[1em]
\left(\mathbf{b} -\mathbf{A}\mathbf{x}^k\right)^T\mathbf{r}^{k} - \alpha^k\left(\mathbf{A}\mathbf{r}^k\right)^T\mathbf{r}^{k}=0\\[1em]
\left(\mathbf{b} -\mathbf{A}\mathbf{x}^k\right)^T\mathbf{r}^{k} = \alpha^k\left(\mathbf{A}\mathbf{r}^k\right)^T\mathbf{r}^{k}
$$

Have a look at the first term on the left-hand side, we have $\mathbf{b} -\mathbf{A}\mathbf{x}^k$. This is nothing else than our search direction, or residual, at iteration $k$, i.e. $\mathbf{r}^k$. Thus, we can simplify this equation as:

$$
\left(\mathbf{r}^k\right)^T\mathbf{r}^{k} = \alpha^k\left(\mathbf{A}\mathbf{r}^k\right)^T\mathbf{r}^{k}
\tag{eq:steepest-decent-alpha-almost-done}
$$

Now let's have a look at the right-hand side. When we have a matrix vector product, i.e. $\mathbf{A}\mathbf{r}^k$, it produces a column vector (for example $\mathbf{c}$) as:

$$
\begin{bmatrix}
a_{11} & a_{12} & \cdots & a_{1,n-1} & a_{1n} \\
a_{21} & a_{22} & \cdots & a_{2,n-1} & a_{2n} \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
a_{n-1,1} & a_{n-1,2} & \cdots & a_{n-1,n-1} & a_{n-1,n} \\
a_{n1} & a_{n2} & \cdots & a_{n,n-1} & a_{nn}
\end{bmatrix}
\begin{bmatrix}
r^{k}_{1} \\
r^{k}_{2} \\
\vdots \\
r^{k}_{n-1} \\
r^{k}_{n}
\end{bmatrix}=
\begin{bmatrix}
c^{k}_{1} \\
c^{k}_{2} \\
\vdots \\
c^{k}_{n-1} \\
c^{k}_{n}
\end{bmatrix}
$$

Then, we have to take the transpose of $\left(\mathbf{A}\mathbf{r}^k\right)^T$, i.e. $\left(\mathbf{c}^k\right)^T$, and multiply that by $\mathbf{r}^k$, i.e.

$$
\left(\mathbf{c}^k\right)^T\mathbf{r}^k\\[1em]
\begin{bmatrix}
c^{k}_{1} & c^{k}_{2} & \dots & c^{k}_{n-1} & c^{k}_{n}
\end{bmatrix}
\begin{bmatrix}
r^{k}_{1} \\
r^{k}_{2} \\
\vdots \\
r^{k}_{n-1} \\
r^{k}_{n}
\end{bmatrix}=
c^{k}_{1}\cdot r^{k}_{1} + c^{k}_{2}\cdot r^{k}_{2} + ... + c^{k}_{1}\cdot r^{k}_{n-1} + c^{k}_{n}\cdot r^{k}_{n}
$$

Again, a row vector multiplied with a column vector is just the scalar product. So, we can also write $\mathbf{c}^k\cdot\mathbf{r}^k=c_i^k\cdot r_i^k$. Since multiplication is commutative (the order of multiplication does not matter), we can write $\mathbf{c}^k\cdot\mathbf{r}^k=\mathbf{r}^k \cdot \mathbf{c}^k$. Writing this again as a row and column vector multiplication, we have:

$$
\mathbf{r}^k \cdot \mathbf{c}^k = \left(\mathbf{r}^k\right)^T \mathbf{c}^k=\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k
$$

We can insert this result back into our equation for $\aplha$, i.e. Eq.(\ref{eq:steepest-decent-alpha-almost-done}), which gives us:

$$
\left(\mathbf{r}^k\right)^T\mathbf{r}^{k} = \alpha^k\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k
$$

Now we divide by $\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k$ and obtain an equation for $\alpha^k$ as:

$$
\alpha^k = \frac{\left(\mathbf{r}^k\right)^T\mathbf{r}^{k}}{\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k}
$$

Even though we have now vectors and matrices in the numerator and denominator, these will all evaluate to scalars, so that $\alpha^k$ will remain a scalar.

You might still be hung up on why we bothered to waste our energy to replace $\left(\mathbf{A}\mathbf{r}^k\right)^T\mathbf{r}^{k}$ with $\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k$. Why did we do that? Was there any point? No, not really, but this is, again, how you will find $\alpha$ introduced in the literature. I am sure who ever came up with this form was just flexing, wanting to show how well they can manipulate linear algebra equations. But there is no need to do so, we could have jsut as well written:

$$
\alpha^k = \frac{\left(\mathbf{r}^k\right)^T\mathbf{r}^{k}}{\left(\mathbf{A}\mathbf{r}^k\right)^T\mathbf{r}^{k}}
$$

We can also get rid of the row vectors now (transposing our column vectors) and write everything as dot products:

$$
\alpha^k = \frac{\mathbf{r}^k\cdot \mathbf{r}^{k}}{\left(\mathbf{A}\mathbf{r}^k\right)\cdot \mathbf{r}^{k}}
$$

Or even:

$$
\alpha^k = \frac{\mathbf{r}^k\cdot \mathbf{r}^{k}}{\mathbf{r}^k\cdot \mathbf{A}\mathbf{r}^k}
$$

All of these equations for $\alpha$ are equivalent and we can use any of them, but, if you open a CFD textbook, you most likely will find $\alpha$ given as:

$$
\alpha^k = \frac{\left(\mathbf{r}^k\right)^T\mathbf{r}^{k}}{\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k}
\tag{eq:steepest-decent-alpha}
$$

This is our final result for alpha. Now that we have $\mathbf{x}^k$, $\mathbf{r}^k$, and $\alpha^k$, we can update our solution using the line search algorithm as:

$$
\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k\mathbf{r}^k
$$

We repeat this until $\mathbf{x}^{k+1}-\mathbf{k}\lt\text{tolerance}$, where we have to specify a suitable tolerance (e.g. $10^{-6}$, something small). Or, we just let it run for a few iterations and *hope* it converges. Ah yes, there it is again, *hope* ...

Thus, we can formulate the Steepest Descent algorithm as follows:

1. Set $\mathbf{x}^k$ from the initial solution or previous timestep/iteration.
2. Compute a search direction $\mathbf{r}^k=\mathbf{b}-\mathbf{Ax}^k$.
3. Compute the scalar multiplier $\alpha^k=(\left(\mathbf{r}^k\right)^T\mathbf{r}^{k})/(\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k)$
4. Find a new solution using the line search algortihm $\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k\mathbf{r}^k$.
5. Go back to step 2 and repeat until $\mathbf{x}^{k+1}-\mathbf{k}\lt\text{tolerance}$, or until the maximum numbers of iterations has been found.

To show you how simple it is to implement the Steepest Descent algorithm, have a look at the following Python implementation using numpy (if you are not familiar with numpy, think of ```@``` as a vector-vector and matrix-vector multiplication operator):

```python
import numpy as np
import matplotlib.pyplot as plt

A = np.array([[2, 1], [1, 2]])
b = np.array([5, 4])

xnew = np.array([-4, -2])
solutions = []
solutions.append(xnew)

for k in range (0, 5):
    xold = xnew
    r = b - A @ xold
    alpha = r.T @ r / (r.T @ A @ r)
    xnew = xold + alpha * r
    solutions.append(xnew)

# print solution table
print('iter |   x    |   y    |')
print('-----|--------|--------|')
for k, x in enumerate(solutions):
    print(f'{k+1:4d} | {x[0]:6.3f} | {x[1]:6.3f} |')
```

The code produces a table at the end to show the values of the new solution we compute. The table is shown below:

```text
iter |   x    |   y    |
-----|--------|--------|
   1 | -4.000 | -2.000 |
   2 |  1.041 |  2.033 |
   3 |  1.905 |  0.953 |
   4 |  1.985 |  1.016 |
   5 |  1.999 |  0.999 |
   6 |  2.000 |  1.000 |
```

The first iteration is our starting point, i.e. let's say that our solution at the previous timestep/iteration was $\mathbf{x}=[-4,-2]$. Then, after only one iteration, we have found a pretty good solution already. I am running this algorithm here for 5 iterations, I am not checking if $\mathbf{x}^{k+1}-\mathbf{x}\lt\text{tolerance}$, because I wanted to show you how quickly this algorithm converges.

Now, think about it for a second. We started with $\mathbf{Ax}=\mathbf{b}$. What we want to do is find a solution for $\mathbf{x}$. The naive approach would be to form the inverse of $\mathbf{A}$ and multiplying both sides by it, which would give us: $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$. But, as we saw, inverting the matrix, using Gaussian elimination, for example, takes too long.

OK, then, in our quest to finding a better algorithm, we looked at LU decomposition, which still uses Gaussian elimination, but stores the result in such a way that we can reuse them. But, even LU decomposition was too slow. Then, we looked at iterative methods, like the Jacobi and the Gauss-Seidel method. While we did not use them to solve a system of linear equations, I said that they are typically very slow to converge, with Gauss-Seidel converging about 2 times faster than Jacobi.

But now, we have arrived at the method of Steepest Descent, where within a single iteration, we get a pretty good estimation of the final solution. After only 5 iterations we have found pretty much the exact solution. And, even if our initial guess for $\mathbf{x}^k$ is far off, this method still converges rapidly. For example, if we multiply our initial solution by a factor of 100 (we have $\mathbf{x}^k=[-400, -200]^T$), then this is the table for 5 iterations (the first iteration being the initial solution):

<!-- wp:kevinbatdorf/code-block-pro {"code":"iter |   x    |   y    |\n\u002d\u002d\u002d\u002d-|\u002d\u002d\u002d\u002d\u002d\u002d\u002d\u002d|\u002d\u002d\u002d\u002d\u002d\u002d\u002d\u002d|\n   1 | -400.0 | -200.0 |\n   2 | -62.25 | 70.20 |\n   3 | -4.355 | -2.177 |\n   4 |  0.984 |  2.094 |\n   5 |  1.900 |  0.950 |\n   6 |  1.984 |  1.017 |","codeHTML":"\u003cpre class=\u0022shiki dark-plus\u0022 style=\u0022background-color: #1E1E1E\u0022 tabindex=\u00220\u0022\u003e\u003ccode\u003e\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003eiter |   x    |   y    |\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e\u002d\u002d\u002d\u002d-|\u002d\u002d\u002d\u002d\u002d\u002d\u002d\u002d|\u002d\u002d\u002d\u002d\u002d\u002d\u002d\u002d|\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e   1 | -400.0 | -200.0 |\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e   2 | -62.25 | 70.20 |\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e   3 | -4.355 | -2.177 |\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e   4 |  0.984 |  2.094 |\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e   5 |  1.900 |  0.950 |\u003c/span\u003e\u003c/span\u003e\n\u003cspan class=\u0022line\u0022\u003e\u003cspan style=\u0022color: #D4D4D4\u0022\u003e   6 |  1.984 |  1.017 |\u003c/span\u003e\u003c/span\u003e\u003c/code\u003e\u003c/pre\u003e","language":"plaintext","theme":"dark-plus","bgColor":"#1E1E1E","textColor":"#D4D4D4","fontSize":".875rem","fontFamily":"Code-Pro-JetBrains-Mono","lineHeight":"1.25rem","clampFonts":false,"lineNumbers":true,"headerType":"none","disablePadding":false,"footerType":"none","enableMaxHeight":false,"seeMoreType":"","seeMoreString":"","seeMoreAfterLine":"","seeMoreTransition":false,"seeMoreCollapse":false,"seeMoreCollapseString":"","highestLineNumber":8,"highlightingHover":false,"lineHighlightColor":"rgba(234, 191, 191, 0.2)","copyButton":true,"copyButtonType":"heroicons","copyButtonUseTextarea":true,"useTabs":false} -->
<div class="wp-block-kevinbatdorf-code-block-pro cbp-has-line-numbers" data-code-block-pro-font-family="Code-Pro-JetBrains-Mono" style="font-size:.875rem;font-family:Code-Pro-JetBrains-Mono,ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;--cbp-line-number-color:#D4D4D4;--cbp-line-number-width:calc(1 * 0.6 * .875rem);line-height:1.25rem;--cbp-tab-width:2;tab-size:var(--cbp-tab-width, 2)"><span role="button" tabindex="0" style="color:#D4D4D4;display:none" aria-label="Copy" class="code-block-pro-copy-button"><pre class="code-block-pro-copy-button-pre" aria-hidden="true"><textarea class="code-block-pro-copy-button-textarea" tabindex="-1" aria-hidden="true" readonly>iter |   x    |   y    |
-----|--------|--------|
   1 | -400.0 | -200.0 |
   2 | -62.25 | 70.20 |
   3 | -4.355 | -2.177 |
   4 |  0.984 |  2.094 |
   5 |  1.900 |  0.950 |
   6 |  1.984 |  1.017 |</textarea></pre><svg xmlns="http://www.w3.org/2000/svg" style="width:24px;height:24px" fill="none" viewBox="0 0 24 24" stroke="currentColor" stroke-width="2"><path class="with-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2m-6 9l2 2 4-4"></path><path class="without-check" stroke-linecap="round" stroke-linejoin="round" d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2"></path></svg></span><pre class="shiki dark-plus" style="background-color: #1E1E1E" tabindex="0"><code><span class="line"><span style="color: #D4D4D4">iter |   x    |   y    |</span></span>
<span class="line"><span style="color: #D4D4D4">-----|--------|--------|</span></span>
<span class="line"><span style="color: #D4D4D4">   1 | -400.0 | -200.0 |</span></span>
<span class="line"><span style="color: #D4D4D4">   2 | -62.25 | 70.20 |</span></span>
<span class="line"><span style="color: #D4D4D4">   3 | -4.355 | -2.177 |</span></span>
<span class="line"><span style="color: #D4D4D4">   4 |  0.984 |  2.094 |</span></span>
<span class="line"><span style="color: #D4D4D4">   5 |  1.900 |  0.950 |</span></span>
<span class="line"><span style="color: #D4D4D4">   6 |  1.984 |  1.017 |</span></span></code></pre></div>
<!-- /wp:kevinbatdorf/code-block-pro -->

While we don't get *perfect* results after 5 iterations, we would likely converge within a few extra iterations. Thus, not only is this method rapid, but it is also pretty insensitive to the initial condition.

But let's go back to our first example, where we started at $\mathbf{x}^k=[-4, -2]^T$. We can visualise the solution vector $\mathbf{x}$ at each iteration and show how the Steepest Descent converges quickly towards the global minimum:

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\steepest_decent_visualised.png" alt="This figure shows the contours of f which we are trying to minimise, and shows the search directions on these contours. The left plot shows a larger portions of the contours of f, while the right plot is a zoomed-in version, shown the algorithm near the minimum of f only." class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

You might be saying, well, didn't we say that the Steepest Descent algorithm takes search directions at a normal angle? Yes, we did say that, but we only used that as a condition to derive $\alpha$. We never actually enforced this condition when we compute a new search direction $\mathbf{r}^{k+1}$. Instead, we simply evaluated the residual at the beginning of each iteration $k$ and set $\mathbf{r}^k=\mathbf{b}-\mathbf{Ax}^k$. That was all.

As a result, we don't get perfectly orthogonal (normal) search directions. As it turns out, the Steepest Descent algorithm does indeed find the local minimum quickly, but, close to the minimum, it tends to oscillate. Once the method is within the vicinity of the exact solution, it struggles to converge. It's like a sprinter who can do the first 90 meters blisteringly fast, but struggles on the last 10 meters in a 100 metre sprint.

If you go back to our example where we started at $\mathbf{x}^k=[-4, -2]^T$, we look at how the solution for the y-coordinate oscillates around the exact solution of 1. It does eventually find it, but there are oscillations. These oscillations arise as we limit our new search direction to be orthogonal to our old search direction. This brute force approach works well if we are far away from the minimum. Close to the minimum, however, things change. Here, we need a different strategy to get quicker convergence.

##### Summary of the Steepest Descent method

In summary, the Steepest Descent method can be written down as follows:

1. Loop from $k=0$ to $N$
2. Compute residual: $\mathbf{r}^{k} = \mathbf{b} -\mathbf{Ax}^{k}$
3. Compute search direction length: $\alpha^k = \frac{\left(\mathbf{r}^k\right)^T\mathbf{r}^{k}}{\left(\mathbf{r}^k\right)^T\mathbf{A}\mathbf{r}^k}$
4. Compute new position on $f$: $\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k\mathbf{r}^k$.
5. Compute residual: $res = \mathbf{x}^{k+1} - \mathbf{x}^{k}$ or $res = \mathbf{b} - \mathbf{Ax}^{k+1}$
6. If $||res||_2 \lt tolerance$ or $k=N$, break, otherwise, go back to step 2.
7. End

#### Issues with the Steepest Descent method

Let's investigate why the Steepest Descent algorithm struggles near the minimum. To do that, I have created two surfaces and visualised the search path of the Steepest Descent algorithm on both. This figure is shown below:

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\steepest_decent_with_deformed_surfaces.png" alt="A sketch shown how the Steepest Descent method performs for surfaces without stretching and with stretching." class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

For both surfaces, I am using a right-hand side vector of $\mathbf{b}=[2,1]^T$, and a starting vector of $\mathbf{x}=[-4,-2]$. However, the coefficient matrix $\mathbf{A}$ is different in both cases. For the left and right surface, I am using the following matrices:

$$
\mathbf{A}_\text{left}=
\begin{bmatrix}
1 & 0\\
0 & 1
\end{bmatrix},\qquad
\mathbf{A}_\text{right}=
\begin{bmatrix}
10 & 7\\
7 & 10
\end{bmatrix}
$$

Let's first try to link the matrices back to what we can see in the figures. The left matrix is the identity matrix. Crucially, it does not have any off-diagonal entries and only the main diagonal contains entries. The surface that we are getting is producing perfectly round iso-contours. What does that mean? Well, at any given point on this surface, if we take the derivative, it will point straight to the minimum.

How about the matrix for the right figure? It does contain off-diagonal elements, and these introduce stretching. These off-diagonal entries stretch our otherwise circular iso-contours to elliptical iso-contours. So, how can we link this now to the performance of the Steepest Descent algorithm? Let's look at the following sketch:

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\circle_ellipse.png" alt="A sketch showing a circle and an ellipse, with tangents drawn on the circle and ellipse. A normal through the tangent will always point to the center of the circle, but a normal through the tangent on an ellipse will not point to the center, as shown in this sketch." class="wp-image-5550" style="width:600px"/></figure>
<!-- /wp:image -->

Here, I show a circle and an ellipse, and I have drawn a tangents at an arbitrary point along the circle and ellipse. If we take the gradient at the intersection of the tangent with the cericle and ellipse, it will point in a direction that is normal to the tangent. We can see that this points straight to the center for the circle, and, it doesn't matter where we were to draw the tangent; the normal to any tangent will always point to the center of the circle!

But how about the ellipse? Well, certainly, there are some places where the normal to a tangent points straight to the center, as shown on the far left of the ellipse. But, these are exceptions, and in most cases, the normal to a tangent (the derivative) points to a location away fromt he center.

As a result, if we use the Steepest Descent method on a surface without any stretching (the matrix $\mathbf{A}$ only contains diagonal entries), then finding the minimum is quick. But, once we have stretching (off-diagonal entries), the Steepest Descent method struggles to point to the center of the surface. And, the stronger the stretching, the worse the performance of the Steepest Descent algorithm gets.

We can also show this mathematically. We looked at it about earlier in this article in Eq.(\ref{eq:condition-number}), repeated below for convenience

$$
\kappa(\mathbf{A})=\frac{\lambda_{max}}{\lambda_{min}}
$$

This depends on the eigenvalues. Since our matrices are 2-by-2, we can quickly evaluate the eigenvalues. For the left matrix, we have:

$$
\det(\mathbf{A}_\text{left}-\lambda\mathbf{I})=
\begin{vmatrix}
1 - \lambda & 0\\
0 & 1 - \lambda
\end{vmatrix}=
(1-\lambda)^2=0
$$

Since we have a 2-by-2 matrix, we expect two eigenvalues. However, there is only one eigenvalue that satisfies $(1-\lambda)^2=0$ and that is $\lambda=1$. Therefore, we have $\lambda_{min}=\lambda_{max}=1$. What does that mean for our condition number? We have:

$$
\kappa(\mathbf{A}_\text{left})=\frac{\lambda_{max}}{\lambda_{min}}=\frac{1}{1}=1
$$

When we discussed the condition number briefly, I mentioned that the closer the eigenvalues are together, the faster the convergence will be. In this case, all eigenvalues (well, all two of them) are exactly the same, and so, we would expect fast convergence. Indeed, looking at the figure above, where we see the surface and the path of the Steepest Descent algorithm, we can see that it converges within a single step.

This is because $\mathbf{A}_\text{left}$ has no off-diagonal entries, and so its contours are perfectly circular, and so no matter where we are on the surface $f$, the gradient will always point in the direction of the center, as we have established by looking at the circle, its tangents and the normals to the tangents.

Now, let's compute the eigenvalues for the second matrix. Here, we have:

$$
\det(\mathbf{A}_\text{right}-\lambda\mathbf{I})=
\begin{vmatrix}
10 - \lambda & 7\\
7 & 10 - \lambda
\end{vmatrix}=
(10-\lambda)^2-7^2=0
$$

From this, we can now compute $\lambda$ to be:

$$
(10-\lambda)^2-7^2=0\\
(10-\lambda)(10-\lambda)-49=0\\
100 - 10\lambda - 10\lambda + \lambda^2 - 49 = 0\\
\lambda^2 -20\lambda + 51 = 0
$$

Using the abc formula, we can find the roots of this expression to be:

$$
\lambda_{1,2}=\frac{-(-20)\pm\sqrt{(-20)^2-4\cdot 1 \cdot 51}}{2\cdot 1} = \frac{20\pm\sqrt{400-204}}{2}=\frac{20\pm\sqrt{196}}{2}=\frac{20\pm 14}{2}\\[1em]
\lambda_1 =\frac{20 + 14}{2} = \frac{34}{2}=17\\[1em]
\lambda_2 =\frac{20 - 14}{2} = \frac{6}{2}=3
$$

Thus, we have $\lambda_{max}=17$ and $\lambda_{min}=3$. What does that do to our condition number?

$$
\kappa(\mathbf{A}_\text{right})=\frac{\lambda_{max}}{\lambda_{min}}=\frac{17}{3}=5.67
$$

Due to the off-diagonal entries, we now get stretching on our surface $f$, and this increases our condition number. A larger condition number means we have to wait longer for the algorithm to converge, and, in the case of the Steepest Descent method, we see that this is because the method has to now zig-zag to the minimum. The stretching means derivatives point no longer to the minimum everywhere on the surface, necessitating this zig-zag approach.

If I borrow the Poisson equation again, which we solve, for example, for an incompressible flow ([SIMPLE](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-solve-incompressible-and-compressible-flows-in-cfd#aioseo-pressure-projection-based-simple), [PISO](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-solve-incompressible-and-compressible-flows-in-cfd#aioseo-how-to-go-from-the-simple-algorithm-to-the-piso-algorithm), [Pressure Projection](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-solve-incompressible-and-compressible-flows-in-cfd#aioseo-pressure-projection-based) methods), then we can see how these matrices with off-diagonal elements (and thus stretched surfaces $f$) come about. The Poisson equation for the pressure may be discretised as:

$$
\frac{1}{\rho}\nabla^2 p^{n+1} = S\\[1em]
\nabla^2 p^{n+1} = \rho S\\[1em]
\frac{\partial^2 p^{n+1}}{\partial x^2} = \rho S\\[1em]
\frac{p_{i+1}^{n+1}-2p_i^{n+1} + p_{i-1}^{n+1}}{(\Delta x)^2} = \rho S_i \\[1em]
p_{i+1}^{n+1}\left[\frac{1}{(\Delta x)^2}\right] + p_{i}^{n+1}\left[\frac{-2}{(\Delta x)^2}\right] + p_{i-1}^{n+1}\left[\frac{1}{(\Delta x)^2}\right] = \rho S_i
p_{i+1}^{n+1} a_E + p_{i}^{n+1} a_P + p_{i-1}^{n+1} a_W = \rho S_i
$$

We can also write this in matrix form as:

$$
\begin{bmatrix}
a_P & a_E & 0 & \cdots & 0 \\
a_W & a_P & a_E & \cdots & 0 \\
0 & a_W & a_P & a_E & 0 \\
\vdots & \vdots & \vdots & \ddots & a_E \\
0 & 0 & 0 & a_W & a_P
\end{bmatrix}
\begin{bmatrix}
p_1 \\
p_2 \\
p_3 \\
\vdots \\
p_n
\end{bmatrix}=
\begin{bmatrix}
\rho S_1 \\
\rho S_2 \\
\rho S_3 \\
\vdots \\
\rho S_n
\end{bmatrix}
$$

If we set $\Delta x = 1$ for simplicity, then we can compute our coefficients as $a_W=a_E=1/(\Delta x)^2=1$ and $a_P=-2/(\Delta x^2)=-2/1=-2$. Our coefficient matrix becomes:

$$
\begin{bmatrix}
-2 & 1 & 0 & \cdots & 0 \\
1 & -2 & 1 & \cdots & 0 \\
0 & 1 & -2 & 1 & 0 \\
\vdots & \vdots & \vdots & \ddots & 1 \\
0 & 0 & 0 & 1 & -2
\end{bmatrix}
$$

It is a symmetric, and positive definite matrix. Crucially, what we can take away from this is that if our stencil of our numerical scheme is symmetric, then our coefficient matrix will be symmetric as well. That is, if our numerical scheme uses $\phi_{i+1}$, we require that we also use $phi_{i-1}$, and the coefficients for both $\phi_{i+1}$ and $\phi_{i-1}$ have to be the same (these will become the entries on the off-diagonals, i.e. $a_E$ and $a_W$).

So, for example, if we used an upwind scheme, it will not use a symmetric stencil, that is, we use either values from $i+1$ *or* $i-1$ in our discretisation, but never both sides. And, since our momentum equations are typically discretised with some form of upwind schemes, we typically have non-symmetric coefficient matrices for our momentum equations and thus can't use the Conjugate Gradient method here.

I should clarify here that even if we don't use the upwind scheme directly, other popular schemes like MUSCL and WENO are all based on the idea of upwinding (i.e. their stencil adjusts based on the direction of the flow).

Returning to our Poisson equation coefficient matrix, if we make this a specific, 2-by-2 matrix again, and compute the eigenvalues for this matrix, we get:

$$
\det(\mathbf{A}_\text{right}-\lambda\mathbf{I})=
\begin{vmatrix}
-2 -\lambda & 1\\
1 & -2 -\lambda
\end{vmatrix}=
(-2 -\lambda)^2-1^2=0
$$

We can again compute $\lambda$ as:

$$
(-2 -\lambda)^2-1^2=0\\
(-2 -\lambda)(-2 -\lambda)-1=0\\
-2 + 2\lambda + 2\lambda + \lambda^2 - 1 = 0\\
\lambda^2 + 4\lambda - 3 = 0
$$

Using the abc formula as before, we can find the roots as:

$$
\lambda_{1,2}=\frac{-(4)\pm\sqrt{(4)^2-4\cdot 1 \cdot (-3)}}{2\cdot 1} = \frac{-4\pm\sqrt{16-(-12)}}{2}=\frac{-4\pm\sqrt{28}}{2}\approx\frac{-4\pm 5.3}{2}\\[1em]
\lambda_1 =\frac{-4 + 5.3}{2} = \frac{1.3}{2}=0.65\\[1em]
\lambda_2 =\frac{-4 - 5.3}{2} = \frac{-9.3}{2}=-4.65
$$

To compute the condition number (which is always positive, that is just its definition), we have to take the absolute values of our eigenvalues, so we get $\lambda_{max}=4.65$ and $\lambda_{min}=0.65$ The condition number is then found to be:

$$
\kappa(\mathbf{A})=\frac{|\lambda_{max}|}{|\lambda_{min}|}=\frac{|-4.65|}{|0.65|}=7.15
$$

So, we can see how a typical example from a pressure Poisson solver, which we frequently use in CFD, especially for incompressible flows, results in matrices that are not ideal for the Steepest Descent method. This is then, where the Conjugate Gradient method shines, so let's have a look at it next.

#### Putting all together: the Conjugate Gradient method at last

Let's approach this problem like an engineer. We need to establish what works really well for the Steepest Descent method, and what doesn't. We retain the favourable properties and try to improve upon the weaknesses.

The Steepest Descent method gets us to the minimum quickly. This is really good, and we would like to use that in our Conjugate Gradient method. However, once we are close to the minimum, it starts oscillating, so the search direction update is not ideal. We want to improve that.

So, we need to consider how we can obtain $\alpha$ in a *better* way. To do that, let's consider the following figure:

<!-- wp:image {"width":"500px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\conjugate_direction.png" alt="A sketch of the surface f and how a possible path on that surface may lead to the minimum. Shown here are orthogonal vectors r to each search direction, and the direct distance to the center of the surface f, indicated by the letter e." class="wp-image-5550" style="width:500px"/></figure>
<!-- /wp:image -->

Here, $\mathbf{x}^k$ a possible location we have obtained on our surface $f$ after trying to get close to the minimum. $\mathbf{x}^{k+1}$ is the location at the next iteration. For both positions $\mathbf{x}^k$ and $\mathbf{x}^{k+1}$, I have also shown the direction vector $\mathbf{r}^k$ and $\mathbf{r}^{k+1}$ in red, which are orthogonal/normal to $\mathbf{x}^k$ and $\mathbf{x}^{k+1}$. This is the direction we would go in the Steepest Descent algorithm.

In blue, I have shown the direct path from $\mathbf{x}^k$ and $\mathbf{x}^{k+1}$ to the center of the surface $f$, indicated by the labels $\mathbf{e}^k$ and $\mathbf{e}^{k+1}$, and it is this direction we ideally would like to take, but we don't know how to get these directions. But, we can express them mathematically as:

$$
\mathbf{e}^k=\mathbf{x}^k - \mathbf{x}^*\\[1em]
\mathbf{e}^{k+1}=\mathbf{x}^{k+1} - \mathbf{x}^*
$$

Here, $\mathbf{x}^*$ is is the location where $f$ has its minimum, as indicated int he figure above as well. That is, the is the location we want to find to get our solution vector $\mathbf{x}$. Therefore, we can also interpret $\mathbf{e}$ as the error between the current location $\mathbf{x}$ and the location where $f$ has its minimum.

As we can see, $\mathbf{x}$ and $\mathbf{e}$ are not necessarily orthogonal to each other. In the Steepest Descent algorithm, we just say that we pick the next search direction to be oprthogonal to the last one, which will guarantee that we will find the minimum eventually, even if it oscillates near the minimum.

In the Conjugate Gradient method, however, we want to try to get closer to the minimum and so we want to move as close as possible in the direction of $\mathbf{e}$ as possible.

Now, let us try to build up some intuition before we go to the next step. We said that if we had an identity matrix, the Steepest Descent method will converge rapidly. It does so no matter where we take the gradient, it will always point directly to the minimum of the surface.

Only when we start to stretch the surface, by introducing off-diagonal terms in our coefficient matrix $\mathbf{A}$, does the Steepest Descent method start to stutter and provide a zig-zag approach to find the solution. So, this leaves us with two possibilities:

1. If we want to continue to use the Steepest Descent algorithm as it is, then we first need to transform our surface $f$ into a space where our surface that we seek to minimuse becomes non-stretched.
2. If we want to leave the surface as it is, we need to change our search direction and it should include our coefficient matrix $\mathbf{A}$, so that the stretching of the off-diagonal terms influences the search direction.

Let's look at both methods in turn. For the first possibility, we need to transform our surface into a non-stretched version. We said that the surface is stretched because of our coefficient matrix $\mathbf{A}$, which has off-diagonal components. So, if we wanted a non-stretched surface $f$, what we need is a matrix that is ideally an identity matrix. And how do we construct an identity matrix from a given coefficient matrix? Well, we have the following definition:

$$
\mathbf{I}=\mathbf{A}^{-1}\mathbf{A}
$$

So, if we have $\mathbf{A}$ given, all we need is its inverse $\mathbf{A}^{-1}$, and we can construct an identity matrix. We can then multiply our surface $f$ by $\mathbf{A}^{-1}$ and have an easy problem to minimise, where the Steepest Descent method will converge rapidly. But, if we have $\mathbf{A}^{-1}$ given, well, then we can solve $\mathbf{Ax}=\mathbf{b}$ explicitly as:

$$
\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}
$$

So, yes, whiule technically this approach would work, it is the nuclear option and it won't help us. We try to avoid having to solve for $\mathbf{A}^{-1}$, and so a transformation of $f$ is off the table.

Let's look at the second approach, where we said that we keep $f$ as it is, but we try to incorporate $\mathbf{A}$ into our search direction.

In the Steepest Descent algorithm, we said that our search directions $\mathbf{r}^k$ and $\mathbf{r}^{k+1}$ should be orthogonal to each other, so that we have:

$$
(\mathbf{r}^k)^T \mathbf{r}^{k+1} = 0
\tag{eq:steepest-decent-search-direction-comparison}
$$

Now we want to insert $\mathbf{A}$ into this expression. As we have seen before, a matrix multiplied by a vector is going to be a vector, so, if we replaced $\mathbf{r}^{k+1}$ by $\mathbf{A}\mathbf{r}^{k+1}$, we still have a vector. Let us insert that into our search direction update and we get:

$$
(\mathbf{r}^k)^T \mathbf{A}\mathbf{r}^{k+1} = 0
\tag{eq:conjugate-gradient-search-direction-comparison}
$$

We said that Eq.(\ref{eq:steepest-decent-search-direction-comparison}) shows orthogonality between $\mathbf{r}^k$ and $\mathbf{r}^{k+1}$. In analogy, we say that Eq.(\ref{eq:conjugate-gradient-search-direction-comparison}) shows A-orthogonality, that is, $(\mathbf{r}^k)^T$ and $\mathbf{A}\mathbf{r}^{k+1}$ are A-orthogonal to each other. Instead of using the word A-orthogonal, we can say that $(\mathbf{r}^k)^T$ and $\mathbf{A}\mathbf{r}^{k+1}$ are *conjugate* to each other. This is were the first part of the name Conjugate Gradient comes from.

Now we do one more step. In Eq.(\ref{eq:conjugate-gradient-search-direction-comparison}), we replace $\mathbf{r}^{k+1}$ by $\mathbf{e}^{k+1}$. We don't know what $\mathbf{e}^{k+1}$, but we can try to eliminate it by construction equations from known quantities that relate to $\mathbf{e}^{k+1}$. We do that because $\mathbf{e}^{k+1}$ points directly to the minimum of $f$, and so, if we can find an equation for $\mathbf{e}^{k+1}$, even just an approximation, we should be better off than with the Steepest Descent algorithm.

Thus, our vectors become:

$$
(\mathbf{r}^k)^T \mathbf{A}\mathbf{e}^{k+1} = 0
$$

Let's look at that visually. We can visualise the difference between $\mathbf{r}^{k+1}$ and $\mathbf{e}^{k+1}$ as follows:

<!-- wp:image {"width":"400px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\cg_search_direction.png" alt="A comparison between r and e close to the minimum of the surface f, where e always points to the minimum of f while r always points away from the minimum, as it is constrained to be orthogonal to the last search direction r" class="wp-image-5550" style="width:400px"/></figure>
<!-- /wp:image -->

We can see the search direction that brought us to $\mathbf{x}^{k+1}$, i.e. the red arrow $\mathbf{r}^k$. Once we reach $\mathbf{x}^{k+1}$, we can see that the next search direction is $\mathbf{r}^{k+1}$ in the Steepest Descent, being orthogonal to $\mathbf{r}^k$. This does not get us to the minimum. But what about $\mathbf{e}^{k+1}$?

I have shown here a circle with its center at $\mathbf{x}^*$ which intersects $\mathbf{x}^{k+1}$. I have also shown a tangent to the circle, and we can see that $\mathbf{e}^{k+1}$ is normal (orthogonal) to the circle's tangent at $\mathbf{x}^{k+1}$. Thus, for any given point $\mathbf{x}^{k+1}$, $\mathbf{e}^{k+1}$ represents the shortest distance to the minimum. Therefore, the expression $(\mathbf{r}^k)^T \mathbf{A}\mathbf{e}^{k+1} = 0$ gives us the shortest direction from $\mathbf{x}^{k+1}$ to $\mathbf{x}^*$.

This is the main difference between the Steepest Descent method and the Conjugate Gradient algorithm. The only thing left for us to do now is to derive an expression for $\mathbf{e}^{k+1}$. Whenever we are tasked with a problem of expressing an unknown property (e.g. $\mathbf{e}^{k+1}$) by known properties, it helps to write out all of the equations that we do know for the unknow property. So let's do that. We already stated that we can express $\mathbf{e}^{k+1}$ as:

$$
\mathbf{e}^{k}=\mathbf{x}^{k} - \mathbf{x}^*\\[1em]
\mathbf{e}^{k+1}=\mathbf{x}^{k+1} - \mathbf{x}^*
$$

We also know from the Steepest Descent algorithm, that we use a line search to update $\alpha$, i.e. Eq.(\ref{eq:line-search}), which was given as:

$$
\mathbf{x}^{k+1}=\mathbf{x}^k + \alpha^k\mathbf{r}^k
$$

We now change $\mathbf{r}$ to $\mathbf{d}$ to indicate that we use an independent search direction in the conjugate gradient method. Thus, our line search becomes:

$$
\mathbf{x}^{k+1}=\mathbf{x}^k + \alpha^k\mathbf{d}^k
$$

If we now subtract $\mathbf{x}^*$ from both sides, we get:

$$
(\mathbf{x}^{k+1} - \mathbf{x}^*)=(\mathbf{x}^k - \mathbf{x}^*) + \alpha^k\mathbf{d}^k\\[1em]
\mathbf{e}^{k+1}=\mathbf{e}^k + \alpha^k\mathbf{d}^k
$$

We also know that once we have found the exact solution $\mathbf{x}^*$, it must satisfy our original problem:

$$
\mathbf{Ax}^* = \mathbf{b}\\[1em]
\mathbf{b} - \mathbf{Ax}^* = 0 
$$

However, if we are still looking for a solution, and $\mathbf{x}^k\ne\mathbf{x}^*$, we know that there will be a residual, that is:

$$
\mathbf{b} - \mathbf{Ax}^k = \mathbf{r}^k 
$$

This residual $\mathbf{r}^k$ was used as the starting direction for the Steepest Descent algorithm. Since the right-hand side vector $\mathbf{b}$ can be expressed as $\mathbf{b} = \mathbf{Ax}^*$, we can insert that into the previous equation and get:

$$
\mathbf{Ax}^* - \mathbf{Ax}^k = \mathbf{r}^k\\[1em]
\mathbf{A}(\mathbf{x}^* - \mathbf{x}^k) = \mathbf{r}^k
$$

Now we can use the definition for $\mathbf{e}^k$ again, which was: $\mathbf{e}^{k}=\mathbf{x}^{k} - \mathbf{x}^*$. Here, the subtraction is the other way around compared to the previous equation. So, we can multiply the previous equation by $-1$ and get:

$$
(-1)\mathbf{A}(\mathbf{x}^* - \mathbf{x}^k) = (-1)\mathbf{r}^k\\[1em]
\mathbf{A}(-\mathbf{x}^* + \mathbf{x}^k) = -\mathbf{r}^k\\[1em]
\mathbf{A}(\mathbf{x}^k -\mathbf{x}^*) = -\mathbf{r}^k
$$

Now that the subtraction is the same as in $\mathbf{e}^{k}=\mathbf{x}^{k} - \mathbf{x}^*$, we can insert and get:

$$
\mathbf{A}\mathbf{e}^{k} = -\mathbf{r}^k
$$

By the same logic, we can also derive an equation for the error at $k+1$. We start with the residual at $k+1$:

$$
\mathbf{b} - \mathbf{Ax}^{k+1} = \mathbf{r}^{k+1} 
$$

Then, we insert $\mathbf{b} = \mathbf{Ax}^*$ for $\mathbf{b}$ and get:

$$
\mathbf{Ax}^* - \mathbf{Ax}^{k+1} = \mathbf{r}^{k+1}\\[1em]
\mathbf{A}(\mathbf{x}^* - \mathbf{x}^{k+1}) = \mathbf{r}^{k+1} 
$$

We multiply again by $-1$:

$$
(-1)\mathbf{A}(\mathbf{x}^* - \mathbf{x}^{k+1}) = (-1)\mathbf{r}^{k+1}\\[1em]
\mathbf{A}(-\mathbf{x}^* + \mathbf{x}^{k+1}) = -\mathbf{r}^{k+1}\\[1em]
\mathbf{A}(\mathbf{x}^{k+1} -\mathbf{x}^*) = -\mathbf{r}^{k+1}
$$

We can insert $\mathbf{e}^{k+1}=\mathbf{x}^{k+1} - \mathbf{x}^*$ and get:

$$
\mathbf{A}\mathbf{e}^{k+1} = -\mathbf{r}^{k+1}
\tag{eq:cg-beta-derivation-1}
$$

With these definition, we can now attempt to find the scalar multiplier $\alpha$ (i.e. the search distance). We start with our conjugate vectors:

$$
(\mathbf{d}^k)^T\mathbf{A}\mathbf{e}^{k+1}=0
\tag{eq:cg-beta-derivation-2}
$$

Here, I have replaced $(\mathbf{r}^k)^T$ by $(\mathbf{d}^k)^T$. I did that because our goal is to find a search direction that is different to the one taken in the Steepest Descent. However, as we will see later, at our first iteration, we have $\mathbf{r} = \mathbf{d}$, that is, in the first iteration both are the same.

However, the Conjugate Gradient method will provide us with a mechanism to independently tune our search direction so that it points closer to the minimum of $f$, and so, we need to distinguish now between the search direction $\mathbf{d}$ and the residual $\mathbf{r}$ (which was used in the Steepest Descent as the search direction as well).

The first step is to eliminate $\mathbf{e}^{k+1}$. We saw that this can be expressed as $\mathbf{e}^{k+1}=\mathbf{e}^k + \alpha^k\mathbf{d}^k$ from our line search approach. So, let us insert that for $\mathbf{e}^{k+1}$. We get:

$$
(\mathbf{d}^k)^T\mathbf{A}(\mathbf{e}^k + \alpha^k\mathbf{d}^k)=0\\[1em]
(\mathbf{d}^k)^T\mathbf{A}\mathbf{e}^k + \alpha^k(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k=0
$$

Since we have already shown that $\mathbf{A}\mathbf{e}^{k} = -\mathbf{r}^k$, we can insert this as well into our equation and get:

$$
-(\mathbf{d}^k)^T\mathbf{r} + \alpha^k(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k=0
$$

Now we divide by $(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k$ and get:

$$
-\frac{(\mathbf{d}^k)^T\mathbf{r}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k} + \alpha^k = 0
$$

Thus, our scalar multiplier $\alpha$ becomes:

$$
\alpha^k = \frac{(\mathbf{d}^k)^T\mathbf{r}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k}
$$

To start the calculation of $\alpha$, as already alluded to before, we set $\mathbf{r}^k=\mathbf{b}-\mathbf{Ax}^k$ and $\mathbf{r}=\mathbf{d}$ at the first iteration, i.e. $k=0$ (or $k=1$ if you must use Fortran, my condolences ...). In this way we can compute $\alpha^{k=0}$.

If we left it at that and set $\mathbf{r}^{k+1}=\mathbf{d}^{k+1}$, well, then we would have created the Steepest Descent algorithm again, but clearly this isn't our goal. So, we need to now look at how we can update the new search direction $\mathbf{d}^{k+1}$ independent of $\mathbf{r}^{k+1}$.

So, let's remind ourselves what equation we already have that we can use to construct $\mathbf{d}^{k+1}$. We will make use of Eq.(\ref{eq:cg-beta-derivation-1}) and Eq.(\ref{eq:cg-beta-derivation-2}) to begin with, which are repated below for convenience:

$$
\mathbf{A}\mathbf{e}^{k+1} = -\mathbf{r}^{k+1}
(\mathbf{d}^k)^T\mathbf{A}\mathbf{e}^{k+1}=0
$$

Inserting $\mathbf{A}\mathbf{e}^{k+1} = -\mathbf{r}^{k+1}$ into $(\mathbf{d}^k)^T\mathbf{A}\mathbf{e}^{k+1}=0$ results in:

$$
(\mathbf{d}^k)^T(-\mathbf{r}^{k+1})=0
$$

We can also write this as:

$$
(-1)\cdot (\mathbf{d}^k)^T \mathbf{r}^{k+1}=0
$$

Now we multiply both sites by $-1$ and get:

$$
(\mathbf{d}^k)^T \mathbf{r}^{k+1}=0
\tag{eq:cg-beta-derivation-3}
$$

So, the search direction at $k$ is orthogonal to the residual $\mathbf{r}$ at the next iteration $k+1$. Let's look at this schematically. The following figure shows how, at $k$, the search direction $\mathbf{d}$ is orthogonal to $\mathbf{r}$ at $k+1$:

<!-- wp:image {"width":"400px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\conjugate_gradient_new_direction.png" alt="A sketch showing how the new search direction d at k+1 points along what was the error e at k+1, which shows that the conjugate gradient method points towards the minimum and thus can achieve faster convergence than the Steepest Descent method." class="wp-image-5550" style="width:400px"/></figure>
<!-- /wp:image -->

The vector $\mathbf{r}^{k+1}$ may have any length, but I am showing it here multiplied by $\alpha$, so that it reaches the local minimum along that search direction. We know that the location reached by $\alpha\mathbf{r}^{k+1}$ is the location where the Steepest Descent method would send us. In the Conjugate Gradient method, however, we go introduce a correction now, that is in the same direction as the old search direction, i.e. $\mathbf{d}^k$, and we multiply that by some unknown scalar multiplier $\beta$.

The goal of $\beta$, as we can see, is to correct the new search direction $\mathbf{d}^{k+1}$ so that it points towards the minimum of the surface $f$. This is shown by the $\mathbf{d}^{k+1}$ pointing in the same direction as the dash-dotted line, shown in blue, which we introduced as the error $\mathbf{e}^{k}$, i.e. the distance between $\mathbf{x}^k$ and the location of the minimum of $f$. If $\mathbf{d}^{k+1}$ would get us all the way to the minimum of $f$, then we would have $\mathbf{d}=\mathbf{e}$.

But, as we are not guaranteed to move all the way to the minimum of $f$, but instead only a fraction of $\mathbf{e}$, they are not the same. But, since we are moving at least in the dirtection of $\mathbf{e}$, we hope to achieve faster convergence compared to the Steepest Descent method, or, at the very least, to avoid oscillations near the minimum.

OK, so far, we have yet again shifted the problem of finding a new search direction $\mathbf{d}^{k+1}$ by expressing it through a new unknown quantity, in this case, $\beta$. Let us write an equation for $\mathbf{d}^{k+1}$ first, and then see how we can obtain $\beta$. Based on the previous figure, we can express the new search direction as a combination of $\mathbf{r}^{k+1}$ and $\beta\mathbf{d}^k$ as:

$$
\mathbf{d}^{k+1} = \mathbf{r}^{k+1} + \beta\mathbf{d}^k
\tag{eq:cg-beta-derivation-4}
$$

Note: I have removed $\alpha$ here from the equation. But that is not an issue, or a simplification. If the length of $\mathbf{r}^{k+1}$ changes, we simply have to change the distance of $\beta\mathbf{d}^k$ so that we are still arriving on a point somewhere on the blue dash-dotted line. Essentially, $\alpha$ is just determining how large the triangle is that is spanned by $\alpha\mathbf{r}^{k+1}$, $\beta\mathbf{d}^k$, and $\mathbf{d}^{k+1}$.

Now let's look at $\mathbf{d}^{k+1}$ in the figure above. You can see that I have drawn a dotted circle around the center of the minimim of $f$. I have also shown the tangent at $\mathbf{x}^k$, and we can see that the new search direction is orthogonal/normal to that tangent.

As we have discussed before, we can say that $\mathbf{d}^{k+1}$ is A-orthogonal, or conjugate, to $\mathbf{d}^k$. We can express that as:

$$
(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^{k+1} = 0
\tag{eq:cg-beta-derivation-5}
$$

Now, let's try to combine that somehow with Eq.(\ref{eq:cg-beta-derivation-4}). In this equation, we know the right-hand side, except for $\beta$. However, we don't know the new search direction $\mathbf{d}^{k+1}$. So, we have two unknowns in this equation. Since we want to determine $\beta$, we need to get rid of the new search direction $\mathbf{d}^{k+1}$.

If we were to multiply Eq.(\ref{eq:cg-beta-derivation-4}) by $(\mathbf{d}^k)^T\mathbf{A}$, then the left-hand side of Eq.(\ref{eq:cg-beta-derivation-4}) would become $(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^{k+1}$. From Eq.(\ref{eq:cg-beta-derivation-5}), we know that this is zero, and so we can effectively eliminate $\mathbf{d}^{k+1}$, and so reduce the number of unknowns in Eq.(\ref{eq:cg-beta-derivation-4}) from 2 to 1.

So, let do that. Let's multiply Eq.(\ref{eq:cg-beta-derivation-4}) by $(\mathbf{d}^k)^T\mathbf{A}$ and we get:

$$
(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^{k+1} = (\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1} + \beta(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k\\[1em]
0 = (\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1} + \beta(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k
$$

We divide this equation now by $(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k$ and get:

$$
0 = \frac{(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k} + \beta\frac{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k} \\[1em]
0 = \frac{(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k} + \beta\frac{1}{1} \\[1em]
0 = \frac{(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k} + \beta \\[1em]
$$

We can solve this for $\beta$ now and get:

$$
\beta = -\frac{(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k}
\tag{eq:beta-with-matrix-multiplication}
$$

Knowing what $\beta is$, we can now determine the new search direction as:

$$
\mathbf{d}^{k+1} = \mathbf{r}^{k+1} + \beta\mathbf{d}^k
$$

This is it, this is the method of the Conjugate Gradient method. I have already talked a few times where the conjugate part is coming from in the name, but what about the gradient part? We saw that we start the Conjugate Gradient method just like the Steepest Descent method, and then correct the search direction. But, at the very first iteration, we set the search direction $\mathbf{d}$ equal to the residual $\mathbf{r}$.

This is coming from the Steepest Descent method, and, since $\mathbf{r}$ is pointing along the gradient, so will $\mathbf{d}$ at the first iteration. So, our starting direction points along the gradient at whichever point $\mathbf{x}$ we are starting at. Hence, since we go initially in the direction of the gradient but then take conjugate, nor orthogonal search directions, we end up with the name Conjugate Gradient.

We could leave things here and move on to the next topic, but, if you look at the Conjugate Gradient method in the literature, you will see a slightly different formulation for $\beta$. In Eq.(\ref{eq:beta-with-matrix-multiplication}), we had:

$$
\beta = -\frac{(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k}
$$

An entirely equivalent expression is:

$$
\beta = \frac{(\mathbf{r}^{k+1})^T \mathbf{r}^{k+1}}{(\mathbf{r}^k)^T \mathbf{r}^k}
\tag{eq:beta-without-matrix-multiplication}
$$

It is just a refomrulated version of $\beta$. Why do we prefer the second verion? Because it does not have any matrix-vector multiplications, only scalar products (or vector vector multiplications). Thus, the second version is much cheaper to evaluate.

Shall I use the same CFD textbook approach and say *it is easy to see that Eq.(\ref{eq:beta-without-matrix-multiplication})follows directly after some manipulation of Eq.(\ref{eq:beta-with-matrix-multiplication})*?

Nah, I'd like to derive my equations in full. For me it is not easy to see that Eq.(\ref{eq:beta-without-matrix-multiplication}) follows from Eq.(\ref{eq:beta-with-matrix-multiplication}) without seeing the derivation.

We start by deriving a useful relation that we will need in a bit. We know that the residual at iteration $k$ and $k+1$ can be written as:

$$
\mathbf{r}^k = \mathbf{b} - \mathbf{A}\mathbf{x}^k\\[1em]
\mathbf{r}^{k+1} = \mathbf{b} - \mathbf{A}\mathbf{x}^{k+1}
$$

We also know, from our line search algorithm, that a new position $\mathbf{x}$ on the surface $f$ at $k+1$, i.e. for the next iteration, can be found as:

$$
\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k \mathbf{d}^k
$$

We now insert $\mathbf{x}^{k+1}$ into our equation for the residual, i.e. $\mathbf{r}^{k+1} = \mathbf{b} - \mathbf{A}\mathbf{x}^{k+1}$ and get:

$$
\mathbf{r}^{k+1}
= \mathbf{b} - \mathbf{A}(\mathbf{x}^k + \alpha^k \mathbf{d}^k)
$$

We can expand the parenthesis and multiply $\mathbf{A}$ by each term. This gives us:

$$
\mathbf{r}^{k+1} = \underbrace{\mathbf{b} - \mathbf{A}\mathbf{x}^k}_{\mathbf{r}^k} - \alpha^k \mathbf{A}\mathbf{d}^k
$$

I have shown that the first two terms combined are nothing else than the residual at iteration $k$, i.e. $\mathbf{r}^k$. Therefore, we can write this also as:

$$
\mathbf{r}^{k+1} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k
$$

We can now rearrange this equation by bringing $\mathbf{r}^k$ onto the other side and get:

$$
\mathbf{r}^{k+1} - \mathbf{r}^k = - \alpha^k \mathbf{A}\mathbf{d}^k
$$

Let's multiply by $-1$ to get rid of the negative sign:

$$
-\mathbf{r}^{k+1} + \mathbf{r}^k = \alpha^k \mathbf{A}\mathbf{d}^k\\[1em]
\mathbf{r}^k - \mathbf{r}^{k+1}  = \alpha^k \mathbf{A}\mathbf{d}^k
$$

As a last step, let's divide by $\alpha^k$:

$$
\frac{1}{\alpha^k}\left(\mathbf{r}^k - \mathbf{r}^{k+1}\right) = \mathbf{A}\mathbf{d}^k
$$

And, just to make it easier for outselves later, let's also swap sides:

$$
\mathbf{A}\mathbf{d}^k = \frac{1}{\alpha^k}\left(\mathbf{r}^k - \mathbf{r}^{k+1}\right)
\tag{eq:Adk-for-beta-derivation}
$$

Good, we will need this result in just a bit. Returning to the equation of $\beta$, we had:

$$
\beta = -\frac{(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}}{(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k}
$$

In the following derivation we will find equivalent expressions for both the numerator $(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}$ and denominator $(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k$.

Let's start with the denominator, i.e. $(\mathbf{d}^k)^T\mathbf{A}\mathbf{d}^k$. We can find an equivalent expression from the scalar multiplier $\alpha^k$. This was given as:

$$
\alpha^k = \frac{(\mathbf{r}^k)^T \mathbf{r}^k}{(\mathbf{d}^k)^T \mathbf{A}\mathbf{d}^k}
$$

Multiplying by $(\mathbf{d}^k)^T \mathbf{A}\mathbf{d}^k$ we get:

$$
\alpha^k(\mathbf{d}^k)^T \mathbf{A}\mathbf{d}^k = \frac{(\mathbf{r}^k)^T \mathbf{r}^k}{1}
$$

Dividing by $\alpha^k$ we get:

$$
(\mathbf{d}^k)^T \mathbf{A}\mathbf{d}^k = \frac{(\mathbf{r}^k)^T \mathbf{r}^k}{\alpha^k}
\tag{eq:denominator}
$$

Let's take that in for a second. The left-hand side shows $(\mathbf{d}^k)^T \mathbf{A}\mathbf{d}^k$, i.e. the denominator for the expression of $\beta$. But, the right-hand side shown an expression that only depends on a single scalar product (or vector-vector multiplication). This will cost less computational effort to evaluate.

Great, now moving on to the numerator $(\mathbf{d}^k)^T\mathbf{A}\mathbf{r}^{k+1}$. Since $\mathbf{A}$ is symmetric, we have $\mathbf{A}^T=\mathbf{A}$. We can use this identity to write the numerator as:

$$
(\mathbf{d}^k)^T\mathbf{A}^T \mathbf{r}^{k+1}
$$

Now we can use the [rules of transposition](https://en.wikibooks.org/wiki/Linear_Algebra/Addition,_Multiplication,_and_Transpose#Transpose). For two generic matrices (of which one can also be a vector as long as the dimensions match) $\mathbf{A}$ and $\mathbf{B}$, we have:

$$
(\mathbf{AB})^T=\mathbf{B}^T\mathbf{A}^T
$$

Thus, we can use this rule for $(\mathbf{d}^k)^T\mathbf{A}^T \mathbf{r}^{k+1}$ and write this as:

$$
(\mathbf{A}\mathbf{d}^k)^T \mathbf{r}^{k+1}
$$

Therefore, the numerator becomes:

$$
(\mathbf{d}^k)^T \mathbf{A}\mathbf{r}^{k+1} = (\mathbf{A}\mathbf{d}^k)^T \mathbf{r}^{k+1}
\tag{eq:dAr-transpose}
$$

Now we use Eq.(\ref{eq:Adk-for-beta-derivation}), which we derived earlier and which was given as:

$$
\mathbf{A}\mathbf{d}^k = \frac{1}{\alpha^k}\left(\mathbf{r}^k - \mathbf{r}^{k+1}\right)
$$

First, we take the transpose of both sides:

$$
\left(\mathbf{A}\mathbf{d}^k\right)^T = \frac{1}{\alpha^k}\left(\mathbf{r}^k - \mathbf{r}^{k+1}\right)^T
$$

We multiply that by $\mathbf{r}^{k+1}$ and get:

$$
\left(\mathbf{A}\mathbf{d}^k\right)^T \mathbf{r}^{k+1} = \frac{1}{\alpha^k}\left(\mathbf{r}^k - \mathbf{r}^{k+1}\right)^T\mathbf{r}^{k+1}
$$

We can see that the left-hand side is equivalent to Eq.(\ref{eq:dAr-transpose})'s right-hand side, so we can substitute:

$$
(\mathbf{d}^k)^T \mathbf{A}\mathbf{r}^{k+1} = \frac{1}{\alpha^k} (\mathbf{r}^k - \mathbf{r}^{k+1})^T \mathbf{r}^{k+1}
$$

We now expand the parenthesis on the right-hand side and get:

$$
(\mathbf{d}^k)^T \mathbf{A}\mathbf{r}^{k+1} = \frac{1}{\alpha^k} \left[ (\mathbf{r}^k)^T \mathbf{r}^{k+1} - (\mathbf{r}^{k+1})^T \mathbf{r}^{k+1} \right]
\tag{eq:beta-derivation-2}
$$

Now we use the identity that the Conjugate Gradient method inherits from the Steepest Descent method. That is, subsequent residuals are orthogonal to each other:

$$
(\mathbf{r}^k)^T \mathbf{r}^{k+1} = 0
$$

Remember, in the Steepest Descent method, we then use $\mathbf{r}^{k+1}$ as the new search direction, whereas in the Conjugate Gradient method, we correct this so that the new search direction is conjugate (taking the coefficient matrix $\mathbf{A}$ into account). For that reason, we have a separate search direction $\mathbf{d}$ from the residual, though, the property that $(\mathbf{r}^k)^T \mathbf{r}^{k+1} = 0$ is still present in the Conjugate Gradient method.

Therefore, we can eliminate $(\mathbf{r}^k)^T \mathbf{r}^{k+1}$ from Eq.(\ref{eq:beta-derivation-2}) and get:

$$
(\mathbf{d}^k)^T \mathbf{A}\mathbf{r}^{k+1} = \frac{1}{\alpha^k} \left[ - (\mathbf{r}^{k+1})^T \mathbf{r}^{k+1} \right]\\[1em]
(\mathbf{d}^k)^T \mathbf{A}\mathbf{r}^{k+1} = -\frac{1}{\alpha^k} (\mathbf{r}^{k+1})^T \mathbf{r}^{k+1}
\tag{eq:beta-denominator}
$$

OK, so, to summarise, from Eq.(\ref{eq:denominator}) and Eq.(\ref{eq:beta-denominator}), we have gotten the following identities:

$$
(\mathbf{d}^k)^T \mathbf{A}\mathbf{d}^k = \frac{(\mathbf{r}^k)^T \mathbf{r}^k}{\alpha^k}\\[1em]
(\mathbf{d}^k)^T \mathbf{A}\mathbf{r}^{k+1} = -\frac{1}{\alpha^k} (\mathbf{r}^{k+1})^T \mathbf{r}^{k+1}
$$

We can now insert these into our definition for $\beta$, i.e.:

$$
\beta = -\frac{(\mathbf{d}^k)^T \mathbf{A}\mathbf{r}^{k+1}}{(\mathbf{d}^k)^T \mathbf{A}\mathbf{d}^k}
$$

Doing so, we get:

$$
\beta = -\frac{-\frac{1}{\alpha^k} (\mathbf{r}^{k+1})^T \mathbf{r}^{k+1}}{\frac{(\mathbf{r}^k)^T \mathbf{r}^k}{\alpha^k}}
$$

The $\alpha^k$ cancels, and the two minus signs become a plus sign, so we can simplify this as:

$$
\beta = \frac{(\mathbf{r}^{k+1})^T \mathbf{r}^{k+1}}{(\mathbf{r}^k)^T \mathbf{r}^k}
$$

Indeed, we have now found an expression for $\beta$ that does not require any matrix-vector multiplications. Horray! Once we have $\beta$, we can compute the search direction $\mathbf{d}$ at $k+1$ which is conjugate to $\mathbf{d}$ at iteration $k$. We have:

$$
\mathbf{d}^{k+1} = \mathbf{r}^{k+1} + \beta\mathbf{d}^k
$$

So, just like in the Steepest Descent method, let us implement the Conjugate Gradient method and see how it performs. We take the same surface $f$, i.e. the same coefficient matrix $\mathbf{A}$ and right-hand side vector $\mathbf{b}$. We also start from the first initial location $\mathbf{x}=[-4, -2]^T$. Here is the implemented algorithm in python again, using numpy:

```python
# conjugate gradient
import numpy as np
import matplotlib.pyplot as plt

A = np.array([[2, 1], [1, 2]])
b = np.array([5, 4])

xnew = np.array([-4, -2])
solutions = []
solutions.append(xnew)

dnew = b - A @ xnew
r = b - A @ xnew
for i in range (0, 5):
    xold = xnew
    dold = dnew
    r = b - A @ xold
    alpha = dold.T @ r / (dold.T @ A @ dold)
    xnew = xold + alpha * dold

    rnew = r - alpha * A @ dold
    beta = rnew.T @ rnew / (r.T @ r)
    dnew = rnew + beta * dold

    solutions.append(xnew)

# print solution table
print('iter |   x    |   y    |')
print('-----|--------|--------|')
for i, x in enumerate(solutions):
    print(f'{i+1:4d} | {x[0]:6.3f} | {x[1]:6.3f} |')
```

On lines 12 and 13, we can see that both the initial residual $\mathbf{r}$ and search direction $\mathbf{d}$ are equivalent at the start of the algorithm. Thus, lines 16-18 are just the Steepest Descent method. On lines 20-22, we do our correction to the search direction, so that ```dnew``` and ```dold```, i.e. the search directions at $k$ and $k+1$ are conjugate to each other.

This algorithm provides the following results:

```text
iter |   x    |   y    |
-----|--------|--------|
   1 | -4.000 | -2.000 |
   2 |  1.041 |  2.033 |
   3 |  2.000 |  1.000 |
   4 |  2.000 |  1.000 |
   5 |  2.000 |  1.000 |
   6 |    nan |    nan |
```

We have found the correct solution to within 3 digits after only 2 iterations! Compare that with the results we obtained from the Steepest Descent method. There, we had the following table:

```text
iter |   x    |   y    |
-----|--------|--------|
   1 | -4.000 | -2.000 |
   2 |  1.041 |  2.033 |
   3 |  1.905 |  0.953 |
   4 |  1.985 |  1.016 |
   5 |  1.999 |  0.999 |
   6 |  2.000 |  1.000 |
```

While the Steepest Descent method also converged to within 3 digits after 6 iteration, it took 3 times longer than the Conjugate Gradient method to get there. And, in comparison to the Conjugate Gradient method, we can clearly see the typical behaviour of the Steepest Descent method, which is that it oscillates near the minimum of $f$.

There is one more thing: After 6 iterations, we get ```nan``` for the Conjugate Gradient method. Now why is that? Refering back to the code, if we run this with debug information available and step through the code, we will find that, after a few iterations, $\beta$ will evaluate to ```nan```.

It does so because it has a division by $(\mathbf{r}^k)^T \mathbf{r}^k$, i.e. the residual. After only a few iterations, we are so close to the minimum of $f$ so that the residual becomes effectively zero (to within machine precision). Therefore, we have division by zero, which sets $\beta$ equal to ```nan```, and then, all values multiplied by $\beta$ also become ```nan```s.

This is yet another indicator that the Conjugate Gradient method converges extremely fast. In reality, we would check if $\mathbf{x}^{k+1}-\mathbf{x}^k \lt tolerance$ and stop the method before we get division by zero, I just wanted to show how quickly we get to the correct solution (minimum of $f$) and how that compares against the Steepest Descent method.

Let me leave you with this analogy, comparing the Steepest Descent method with the Conjugate Gradient algorithm.

I don't play golf (I suppose for the same reason as [Franky Boyle](https://youtu.be/qft-GHY8w94?si=B6zJvZiFrTjHMU8O&t=785)), but the golf ball is one of those things we can't avoid in fluid mechanics. It is often used to explain the favourable properties of turbulent boundary layers and how they can actually decrease drag (which, the [Mythbusters](https://www.youtube.com/watch?v=VUiGhyHC-1A) put to the extreme to build a dimpled, golf ball-like car to reduce drag), as they carry more energy than laminar boundary layers.

Be that as it may, even I know that a golfer has different clubs available at their disposal. You would use a driver to get the golf ball as close to the green (where the hole is) as possible. The driver gives you a good distance, but not good precision. Once you are on the green, you want to use a putter, which doesn't give you distance, but precision instead.

The Steepest Descent method is like a driver; you can be as far away from the function's minimum (i.e. the green/hole in this analogy) as you want, the Steepest Descent algorithm will get you very quickly into the vicinity of the minimum (i.e. the green). The Conjugate Gradient algorithm uses the same approach as the Steepest Descent method initially for exactly that reason; fast convergence *towards* the minimum of $f$.

Once our golf ball is on the green, a golfer would now switch to a putter, which trades precision for distance. We won't get the golf ball far, but we have a lot more control over its direction of travel.

In our analogy, the Steepest Descent method is like continuing to use a driver, even once we are on the green. The Conjugate Gradient method is like a putter, as we continuisly correct our direction of travel, and so, we have more control (precision) on getting the ball into the hole (finding the minimum of $f$), and therefore converge faster (having a lower handicap).

It's good to see that all of that Nintendo Wii sports resort golf training wasn't in vain ...

OK, thus far, we have made the restriction that all of our matrices ought to be symemtric, but that really means all we can do is solve pressure Poisson equations, or heat diffusion type equations with the Conjugate Gradient method. We can't use it for the momentum equation (unless we employ a symmetric stencil for our non-linear term, and there really is only one scheme that comes to mind that provides us with this property; the [JST scheme](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/space-and-time-integration-schemes-for-cfd-applications#aioseo-jst-scheme)). What if we wanted to generalise the Conjugate Gradient algorithm?

Well, fear not, others have done exactly that for us so that we can use the Conjugate Gradient method for non symmetric matrices as well, and by far the most popular variant (there are several), is the Bi-Conjugate Gradient (Stabilised) version, or BiCGStab/BCGS, depending on which literature/software packages you use. We will review this method in a second as well. But first, let's summarise our Conjugate Gradient method, so that we have an easy cheat sheet should we need to implement it in the future.

##### Summary of the Conjugate Gradient method

In summary, the Conjugate Gradient method can be written down as follows:

1. Loop from $k=0$ to $N$
2. Compute residual: $\mathbf{r}^{k} = \mathbf{b} -\mathbf{Ax}^{k}$
3. Compute search direction length: $\alpha^k = \frac{\left(\mathbf{d}^k\right)^T\mathbf{r}^{k}}{\left(\mathbf{d}^k\right)^T\mathbf{A}\mathbf{d}^k}$
4. Compute new position on $f$: $\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k\mathbf{d}^k$.
5. Compute new residual: $\mathbf{r}^{k+1} = \mathbf{b} -\mathbf{Ax}^{k+1}$
6. Compute correction along $\mathbf{d}$: $\beta^k = \frac{\left(\mathbf{r}^{k+1}\right)^T\mathbf{r}^{k+1}}{\left(\mathbf{r}^k\right)^T\mathbf{r}^k}$
7. Compute new (conjugate) search direction: $\mathbf{d}^{k+1} = \mathbf{r}^{k+1} + \beta^k\mathbf{d}^{k}$
8. Compute residual: $res = \mathbf{x}^{k+1} - \mathbf{x}^{k}$ or $res = \mathbf{b} - \mathbf{Ax}^{k+1}$
9. If $||res||_2 \lt tolerance$ or $k=N$, break, otherwise, go back to step 2.
10. End

### Interlude: The Conjugate Gradient algorithm and the Krylov subspace

Remember our friend Krylov? We spent a lot of time looking at what the Krylov subspace is and how it is constructed. It is time to return to Krylov and look at why the Conjuage Gradient method (along with the Steepest Descent method and those that we are yet to review in this section of the article) are said to be Krylov subspace methods.

To do that, let me take you on a journey. *On occasion*, I have mentioned that I am very much German, and that is true, I was born and raised in there. However, having both German and Norwegian parents, I am *technically speaking* half-German and half-Norwegian (*noen andre nordmenn her?*).

In any case, as a result, I have spend lots of my summer holidays in Norway, visiting family, and one place in particular: Narvik. It's a small little city above the northern polar circle which makes it perfect to see northern lights in winter (which I was lucky enough to see a few times). In summer, though, the sun never sets, and you can enjoy a midnight swim in the Fjord. Here is an arial view onto Narvik:

<!-- wp:image {"width":"600px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\narvik_linken.webp" alt="An arial view onto the city of Narvik, showing the town to the left and the mountains to the right." class="wp-image-5550" style="width:600px"/><figcaption class="wp-element-caption">Figure reproduced from <a href="https://lkab.com/en/community-development/lkab-in-narvik/" target="_blank" rel="noopener" title="">LKAB</a>.</figcaption></figure>
<!-- /wp:image -->

You see the town at the base of the mountain on the left. Ignore the airport on the bottom left, that has long closed now (much to my anoyance, I never got a [Saba](https://en.wikipedia.org/wiki/Juancho_E._Yrausquin_Airport)-style landing in my life ...). The mountain on the right is what we are interested in.

Whenever I visited Narvik as a child/teenager, I would scale the mountain at least once, but often times more than that. My uncle, who is very much into Mountaineering, took me on a 1 week tour through Norway, trekking and climbing some of the larger peaks in the south.

This trip culminated in us sleeing in a hut that was officially ranked Norway's heighest building (build ontop of a mountain of about 2200 meters). They had bunkbeds and I was the only one sleeping on the top, so, for at least a night, I was Norway's heigest person (pun intended).

So, even though climbing mountains was not one of my active hobbys (that was a pretty difficult hobby to exercise to someone living in northern Germany where the heighest peaks around would be some [landfills full our garbage](https://de.wikipedia.org/wiki/M%C3%BCllberg_Hummelsb%C3%BCttel) and recultivated to look pretty), it is fair to say that I had a general interest in it.

So, imagine my delight when we went on a trip to Leines one day and I was gretted with a mountain that was yet to be climbed by me; the Kraktindan:

<!-- wp:image {"width":"800px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\kraktinden.webp" alt="Arial view from the Kraktindan mountain at Leines beach" class="wp-image-5550" style="width:800px"/><figcaption class="wp-element-caption">Figure reproduced from <a href="https://www.alltrails.com/de/route/norway/nordland/topptur-til-kraktindan-pa-leines" target="_blank" rel="noopener" title="">AllTrails</a>.</figcaption></figure>
<!-- /wp:image -->

We actually came here quite a few times, and I was always looking at the mountain, wanting to climb it at least once. But, everytime we came, it was raining, or it rained the day before and the trail was no good. But then, one day, a golden opportunity presented itself. The weather was cloudy, but no rain was forecasted. The trail was in good shape, and so I ascended to the peak.

I didn't know back then that this mountain in particular had a bit of a reputation. No one warned me about it because, well, I didn't really think to tell anyone that I was going there. I think I left a message before leaving, and off I went. The ascend was not too bad, and I must have reached the peak within two or three hours. It wasn't a massive mountain, and as long as you could see the peak, you knew where to go.

But then the mountain showed me why it had developed a reputation amongst locals, and why you should, perhaps, go a bit better prepared, like, havig a map, or even better, a GPS with you (we are talking about the early 2000s here, smartphones did not exist yet and the best part about my mobile at that time was [Snake 2](https://helpfulsheep.com/snake/)).

While it was cloudy for the whole duration of the ascend, once I reached to peak, I realised that I was very close to the clouds. I took a bit of a break, taking in the view, and, without much of a warning, the cloud cover descended, and all of a sudden there was a 5 meter visibility at best. While some footpaths existed near the base of the mountain, there were no clear paths near the peak.

So, here I was, no map, no GPS, no compass, only me and my teenage ignorance for Norwegian climate. I decided, perhaps it was time to go back home, but where to go? I knew roughly which way to go, but I didn't know the exact bearing I should take. On the ascend, I had to go through a dense forrest through which a clear footpath existed, but the last 30-40% of the climb was just walking up the mountain.

So, I needed to go back to the footpath that led me through the forrest, if I could not reach it, I would have no idea how to get through the forrest. I was looking for a small opening in the forrest without knowing where to go, exactly.

So, imagine this is you, how would you proceed? Here is what I did: I didn't have any external elements to rely on, like maps or GPS, so I had to entirely rely on what I could see, feel, or sense. Now, there is a big difference between going up or down a mountain, and so I made sure that with each step that I took, that I would go downhill.

Step-by-step I descended, slowely but steadily. My reasoning was that I would eventually reach the forrest this way, and I could just walk along the bundary of the forrest until I would find the footpath. I managed to do just that and eventually found my way back down back to safety.

So, how does this relate to Krylov? I didn't know it back then, but I was applying a Krylov subspace to my problem to get back to safety.

In the Steepest Descent method and the Conjugate Gradient algorithm, we tried to find the minimum of a surface $f$. I was standing very much on a surface $f$ and I was trying to reach the minimum of $f$ as well (i.e. sealevel). Just like the Steepest Descent method, or the Conjugate Gradient algorithm, I was looking for a path downhill.

So, even though I was at a location $\mathbf{x}^k$ on a surface $f$ (the mountain), I could have taken a step in any direction, i.e. I could have taken any of the $360^\circ$ directions available to me, yet I restricted myself to a subset of directons: any direction that was moving down the mountain.

Now, my subset was not random, they all obeyed similar properties (all directions pointed down the hill). So, not only were the directions I could take a subset of the $360^\circ$, but they also formed a subspace, since all directions had the same physical properties of getting me down the mountain. None of the direction I took to get to $\mathbf{x}^{k+1}$ were taking me uphill, for example.

The same is true with the Steepest Descent method, or the Conjugate Gradient algorithm. With each iteration $k$, we try to descend on $f$ until we have found its minimum. We are not taking any direction that take us away from the minimum, we are always descending. For that reason, the direction we take are taken from a subspace of possible directions that all get us closer to the minimum of $f$.

In a nutshell, this is why we say these methods are part of a Krylov subspace. Now, let's revist the *Tom is always right* algorithm. Remember that abomination? Let me plot it below again:

<!-- wp:image {"width":"500px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\line_search_algorithm.png" alt="This figure shows the same function f as seen previously with contours indicating the levels of f in a range of -5 to 5 for both x and y, with x=2 and y=1 being the location of the minimum of f. A random walk shows how one may get to the minimim of f, indicated by points joined together through red dashed lines." class="wp-image-5550" style="width:500px"/></figure>
<!-- /wp:image -->

Look closely at the path we take here. If we are starting at the point shown on the bottom left, then between point 5 and 6, we are going uphill, all of a sudden. We are not monotonically descending on $f$, and so, if we tolerate this behaviour, it means that this method is not guaranteed to converge. In fact, there is a very real possibility that it will diverge.

So, is this method a Krylov subspace method? Well, yes and no. It is, because the subspace we are now defining is that any direction we take has to point to the right compared to the previous direction. So, we are consistently taking a new direction from this subspace. As we have seen, it is a pretty bad subspace, as the method has the potential to diverge (or even just oscillate and never find the minimum).

However, it is not a Krylov subspace method in the sense that we have looked at so far, i.e. our subspace is one where we take a direction that is guaranteed to move down the surface of $f$, never up. So, in that sense, *Tom is always right* fails, and it is not a Krylov subspace method, at leat in the context we use in CFD.

We can even show this mathematically, if we wanted. Remember the definition of the Krylov subspace? It was given as:

$$
\mathcal{K}(\mathbf{A},\mathbf{v})=\text{span}\{\mathbf{v},\mathbf{Av},\mathbf{A}^2\mathbf{v},\mathbf{A}^3\mathbf{v},...,\mathbf{A}^{r-1}\mathbf{v}\}
$$

Let's use the Steepest Descent method to show it is a Krylov subspace method. We saw that the Conjugate Gradient algorithm takes the Steepest Descent method as a starting point and then corrects it by taking conjugate directions. Thus, if we can show that the Steepest Descent method satisfies this definition, then we can say that the Conjugate Gradient algorithm will too.

The first step was to compute the intiial residual, which we take as our starting direction in both the Steepest Descent method, as well as the Conjugate Gradient algorithm. This was given as:

$$
\mathbf{r}^{k} = \mathbf{b} -\mathbf{Ax}^{k}
$$

At the next iteration, the residual is computed as:

$$
\mathbf{r}^{k+1} = \mathbf{b} -\mathbf{Ax}^{k+1}
$$

We also have our line search update, i.e. an equation that allows us to move from $\mathbf{x}^k$ to $\mathbf{x}^{k+1}$. This was given as:

$$
\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k\mathbf{r}^k
$$

We can insert $\mathbf{x}^{k+1}$ into the equation for $\mathbf{r}^{k+1}$, which gives us:

$$
\mathbf{r}^{k+1} = \mathbf{b} -\mathbf{A}(\mathbf{x}^k + \alpha^k\mathbf{r}^k)\\[1em]
\mathbf{r}^{k+1} = \mathbf{b} -\mathbf{A}\mathbf{x}^k - \alpha^k\mathbf{A}\mathbf{r}^k\\[1em]
\mathbf{r}^{k+1} = \mathbf{r}^{k} - \alpha^k\mathbf{A}\mathbf{r}^k
\tag{eq:steepest-descent-krylov-update}
$$

Using the equality that a matrix raised to the power of zero is just an identity matrix, i.e. $\mathbf{A}^0 = \mathbf{I}$, we could write this also as:

$$
\mathbf{r}^{k+1} = \mathbf{A}^0\mathbf{r}^{k} - \alpha^k\mathbf{A}^1\mathbf{r}^k
$$

So, the first term only depends on $\mathbf{r}^{k}$, while the second term depends on $\mathbf{A}$ and $\mathbf{r}^k$. We could say that $\mathbf{r}^{k+1} = f(\mathbf{A}, \mathbf{r}^k)$, or, we could say that we have formed a Krylov subspace as:

$$
\mathcal{K}(\mathbf{A},\mathbf{r}^k)=\text{span}\{\mathbf{r}^k,\mathbf{Ar}^k\}
$$

Notice how the subspace here is defined for the initial residual $\mathbf{r}^k$, this is key. Each subsequent residual $\mathbf{r}^{k+1},\mathbf{r}^{k+2},\mathbf{r}^{k+3},...,\mathbf{r}^{N}$ will have a dependence on $\mathbf{r}^k$. The Krylov subspace is formed using this initial residual $\mathbf{r}^k$ and so we need to express how $\mathbf{r}^{k+1},\mathbf{r}^{k+2},\mathbf{r}^{k+3},...,\mathbf{r}^{N}$ depend on $\mathbf{r}^k$.

Let's do that for $\mathbf{r}^{k+2}$. First, let's write down the definition:

$$
\mathbf{r}^{k+2} = \mathbf{r}^{k+1} - \alpha^k\mathbf{A}\mathbf{r}^{k+1}
$$

This is just Eq.(\ref{eq:steepest-descent-krylov-update}) written now for $k+2$ instead of $k+1$. We have $\mathbf{r}^{k+1}$ on the right-hand side, but we want to express that in terms of $\mathbf{r}^{k}$. We know how to express $\mathbf{r}^{k+1}$ using only $\mathbf{r}^{k}$ from Eq.(\ref{eq:steepest-descent-krylov-update}). So, let's insert:

$$
\mathbf{r}^{k+2} = \left(\mathbf{r}^{k} - \alpha^k\mathbf{A}\mathbf{r}^k\right) - \alpha^{k+1}\mathbf{A}\left(\mathbf{r}^{k} - \alpha^k\mathbf{A}\mathbf{r}^k\right)
$$

Now we just remove the brackets:

$$
\mathbf{r}^{k+2} = \mathbf{r}^{k} - \alpha^k\mathbf{A}\mathbf{r}^k - \alpha^{k+1}\mathbf{A}\mathbf{r}^{k} + \alpha^k\alpha^{k+1}\mathbf{A}^2\mathbf{r}^k
$$

We can now group terms together and arrive at:

$$
\mathbf{r}^{k+2} = \mathbf{r}^{k} - (\alpha^k + \alpha^{k+1})\mathbf{A}\mathbf{r}^k + \alpha^k\alpha^{k+1}\mathbf{A}^2\mathbf{r}^k
\tag{eq:steepest-descent-krylov-update-2}
$$

Thus, $\mathbf{r}^{k+2}$ depends on $\mathbf{r}^{k}$, $\mathbf{A}\mathbf{r}^k$, and $\mathbf{A}^2\mathbf{r}^k$. We can also write this as:

$$
\mathcal{K}(\mathbf{A},\mathbf{r}^k)=\text{span}\{\mathbf{r}^k,\mathbf{Ar}^k,\mathbf{A}^2\mathbf{r}^k\}
$$

We could do the same now for $\mathbf{r}^{k+3}$ and get:

$$
\mathbf{r}^{k+3} = \mathbf{r}^{k+2} - \alpha^k\mathbf{A}\mathbf{r}^{k+2}
$$

With the definition for $\mathbf{r}^{k+2}$ from Eq.(\ref{eq:steepest-descent-krylov-update-2}), we could insert this and get:

$$
\mathbf{r}^{k+3} = \left(\mathbf{r}^{k} - (\alpha^k + \alpha^{k+1})\mathbf{A}\mathbf{r}^k + \alpha^k\alpha^{k+1}\mathbf{A}^2\mathbf{r}^k\right) - \alpha^k\mathbf{A}\left(\mathbf{r}^{k} - (\alpha^k + \alpha^{k+1})\mathbf{A}\mathbf{r}^k + \alpha^k\alpha^{k+1}\mathbf{A}^2\mathbf{r}^k\right)
$$

Removing the brackets again by multiplying out terms we get:

$$
\mathbf{r}^{k+3} = \mathbf{r}^{k} - (\alpha^k + \alpha^{k+1})\mathbf{A}\mathbf{r}^k + \alpha^k\alpha^{k+1}\mathbf{A}^2\mathbf{r}^k - \alpha^k\mathbf{A}\mathbf{r}^{k} + \alpha^k(\alpha^k + \alpha^{k+1})\mathbf{A}^2\mathbf{r}^k - (\alpha^k)^2\alpha^{k+1}\mathbf{A}^3\mathbf{r}^k
$$

Grouping terms to clean up this equation a abit we get:

$$
\mathbf{r}^{k+3} = \mathbf{r}^{k} - (2\alpha^k + \alpha^{k+1})\mathbf{A}\mathbf{r}^k + \alpha^k(\alpha^k + 2\alpha^{k+1})\mathbf{A}^2\mathbf{r}^k - (\alpha^k)^2\alpha^{k+1}\mathbf{A}^3\mathbf{r}^k
$$

As we can see from the equation of $\mathbf{r}^{k+3}$, we have a Krylov subspace forming as:

$$
\mathcal{K}(\mathbf{A},\mathbf{r}^k)=\text{span}\{\mathbf{r}^k,\mathbf{Ar}^k,\mathbf{A}^2\mathbf{r}^k,\mathbf{A}^3\mathbf{r}^k\}
$$

We can repeat this now for all $N$ iterations as in $k=1,2,3,...,N$, but I don't have time to derive equations to infinity at the moment, so you'll have to trust me that, in general, we can write that the Steepest Descent method forms the following Krylov subspace:

$$
\mathcal{K}(\mathbf{A},\mathbf{r}^k)=\text{span}\{\mathbf{r}^k,\mathbf{Ar}^k,\mathbf{A}^2\mathbf{r}^k,\mathbf{A}^3\mathbf{r}^k,...,\mathbf{A}^N\mathbf{r}^k\}
$$

Bring this back to our mountain analogy; we seek to reach the base of the mountain (minimise $f$) and, as long as I go downhill, I will eventually reach the base. So, at every position $\mathbf{x}^k$ on the mountain (our surface $f$), I have a subspace of directions I can take, all of which are pointing downhill. There are some directions that will point up the mountain (up the surface $f$), and these are not valid directions for us to take.

With each step down the mountain (with each iteration $k=1,2,3,...,N$), we take a path down the hill. The Steepest Descent or Conjugate Gradient method just tell us how to take those steps (i.e. they will tell us the specific direction to take from all the allowable directions within our subspace, i.e. all those directions pointing downhill).

Why is the Krylov subspace dependent on the initial location (or rather, residual), at the first iteration $k$? Well, think about it this way, if you go down a mountain from two different locations, you wouldn't take the same path to go down. There may be overlaps, but each path you taker is unique.

This is what $\mathbf{r}^k,\mathbf{Ar}^k,\mathbf{A}^2\mathbf{r}^k,\mathbf{A}^3\mathbf{r}^k,...,\mathbf{A}^N\mathbf{r}^k$ expresses, i.e. these terms were derived by finding $\mathbf{r}^{k+1}$, $\mathbf{r}^{k+2}$, $\mathbf{r}^{k+3}$, and so on. Since $\mathbf{r}^k$ also depends on the location $\mathbf{x}^k$ as in $\mathbf{r}^{k} = \mathbf{b} -\mathbf{Ax}^{k}$, we can see how each of these terms in the Krylov subspace definitions simply show us the path we take to get to the base of the mountain (minimum of $f$).

Hopefully I was able to make that connection between the Conjugate Gradient algorithm (as well as the Steepest Descent method) and the Krylov subspace. As they say, its not rocket science, just a case of no one (at least in the landscape of CFD textbooks to the best of my knowledge) making an effort to show this link. I know, I spend quite a bit on this method, but I figured it is important. Finally, let us move on to the next method.

### The Bi-Conjugate Gradient method

I have alluded to the fact that the Conjugate Gradient algorithm can only be used for symmetric matrices, which typically limits its application area to [elliptic equations](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/what-are-hyperbolic-parabolic-and-elliptic-equations-in-cfd#aioseo-elliptic-flows). In elliptic equations, the flow has no preferred direction of travel. Pressure waves, for example, travel in all directione equally. For that reason, I can stand anywhere in the same room as you and still hear you talk, regardless of my position inside the room.

Or, we may sit around a fire, it doesn't matter where we sit around the fire, as long as we sit the same distance away from it, we will receive the same amount of heat. Thus, for these applications where the direction of propagation is not important, numerical schemes with symmetrical stencils are used that lead to symmetrical coefficient matrices. We can then use the Conjugate Gradient algorithm to solve such a system of elliptical equations.

The Navier-Stokes equations, however, are the biggest party-poopers known to mankind, and so, as always, we have to find a way to modify our otherwise perfectly fine Conjugate Gradient algorithm so that we can deal with the properties of the Navier-Stokes equations, namely their [hyperbolic behaviour](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/what-are-hyperbolic-parabolic-and-elliptic-equations-in-cfd#aioseo-hyperbolic-flows).

In hyperbolic flows, we have a preferred direction of information propagation. We can compute these directions, which turn out to be our [characteristic directions](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/what-are-hyperbolic-parabolic-and-elliptic-equations-in-cfd#aioseo-precurser) (usually abbreviated to just characteristics). For example, a shock wave forming around an object travelling at supersonic speeds.

These type of equations typically require the direction of travel to be baked into the discretisation scheme, which we collectively group under the term *upwinding*, as there are more than just the [upwind scheme](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/space-and-time-integration-schemes-for-cfd-applications#aioseo-first-order-upwind-scheme) that make use of upwinding itself as a concept.

Take the implicit advecction equation in 1D, for example. We can write this as:

$$
\frac{\phi^{n+1} - \phi^k}{\Delta t} + a\frac{\partial \phi^{n+1}}{\partial x}
$$

If we assume $a>0$, and using a first-order upwind scheme here with a finite difference discretisation, we get:

$$
\frac{\phi^{n+1}_i - \phi^k_i}{\Delta t} + a\frac{\phi^{n+1}_i - \phi^{n+1}_{i-1}}{\Delta x}
$$

Writing this in matrix coefficient form, we have:

$$
\phi^{n+1}_i\left[\underbrace{\frac{1}{\Delta t} + \frac{1}{\Delta x}}_{a_P}\right] + \phi^{n+1}_{i-1}\left[\underbrace{-\frac{1}{\Delta x}}_{a_W}\right] = \frac{\phi^n_i}{\Delta t}
$$

We can also write this in matrix form as:

$$
\mathbf{Ax}=\mathbf{b}\\[1em]
\underbrace{\begin{bmatrix}
a_P & 0 & 0 & \cdots & 0 \\[1em]
a_W & a_P & 0 & \cdots & 0 \\[1em]
0 & a_W & a_P & \ddots & \vdots \\[1em]
\vdots & \ddots & \ddots & \ddots & 0 \\[1em]
0 & \cdots & 0 & a_W & a_P
\end{bmatrix}}_\mathbf{A}
\underbrace{\begin{bmatrix}
\phi_1 \\[1em]
\phi_2 \\[1em]
\phi_3 \\[1em]
\vdots \\[1em]
\phi_{nx}
\end{bmatrix}}_\mathbf{x}=
\underbrace{\begin{bmatrix}
\frac{\phi_1^n}{\Delta t} \\[1em]
\frac{\phi_2^n}{\Delta t} \\[1em]
\frac{\phi_3^n}{\Delta t} \\[1em]
\vdots \\[1em]
\frac{\phi_{nx}^n}{\Delta t}
\end{bmatrix}}_\mathbf{b}
$$

As we can see, the coefficient matrix $\mathbf{A}$ is no longer symmetric, and this is a very representative example for hyperbolic flows, and thus, the Navier-Stokes equations in general. So, if we want to use the Conjugate Gradient method here, we need to modify it to allow for non-symmetric matrices as well.

Now, we could go through all of the derivation again and assume that our coefficient matrix is no longer symmetric. If we did that, we would find that we wouldn't be quite able to find an algorithm, as the assumption of using a symmetric matrix helped us a lot to find an algorithm in the first place. Thus, what researchers have done instead is to find ways to *trick* the Conjugate Gradient method into working with non-symmetric matrices.

Think of the following analogy. You want to go to work and sit down in your car and try to start it. But, your car does not start. It turns out the battery is dead. What do you do? Well, you have two options:

- You can go back home, order a new battery online, and replace the flat battery in your car.
- Or, you could ask one of your neighbours to help you jump start the car vy connecting your battery to their car battery.

In this analogy, buying a new battery is like taking the Conjugate Gradient algorithm and rederiving it assuming we have a non-symmetric matrix (i.e. we fix the part that is broken within the Conjugate Gradient algorithm). Jump starting the car (leaving the broken battery in place) is like ignoring the fact that we can't handle non-symmetric matrices and we try to find a workaround instead (we connect our car's battery to another one).

In the world of CFD, there are two methods that are predominantly used *to jump start the Conjugate Gradient algorithm*. The first is the Bi-Conjugate Gradient method (and its stabilised form, i.e. Bi-Conjugate Gradient Stabilised), and the second is the Generalised Minimal Residual (GMRES). We will look at both, but in this section, we will derive the main idea behind the Bi-Conjugate Gradient method first.

So, let's think about this a bit. We want to fix the Conjugate Gradient algorithm and allow non-symmetric matrices to work, but we know that it can only work with symmetric matrices, what can we do? Well, we know that we can't apply it directly to our non-symmetric matrix, but, perhaps we could modify the system in such a way that it contains the non-symmetric matrix in a way that looks symmetric.

The key idea here is that we embed the non-symmetric matrix $\mathbf{A}$ in a block matrix $\mathbf{H}$. Let's write out a generic block matrix first to see how it is defined:

$$
\mathbf{H} =
\begin{bmatrix}
\mathbf{B} & \mathbf{C} \\
\mathbf{D} & \mathbf{E}
\end{bmatrix}
$$

Here, $\mathbf{B}$, $\mathbf{C}$, $\mathbf{D}$, and $\mathbf{E}$ are all matrices themselves. If we could set one or more matrices equal to $\mathbf{A}$, such that when we transpose the system we get exactly the same system back, then $\mathbf{H}$ is symmetric. The idea here is that we then apply the Conjugate Gradient method to $\mathbf{H}$, instead of $\mathbf{A}$, which works without modifications. So, let's find an arrangement in which we construct a symmetrich matrix $\mathbf{H}$ which only depends on $\mathbf{A}$.

First, let's see how block matrices are tranposed, which is shown in the following:

$$
\mathbf{H}^T =
\begin{bmatrix}
\mathbf{B} & \mathbf{C} \\
\mathbf{D} & \mathbf{E}
\end{bmatrix}^T=
\begin{bmatrix}
\mathbf{B}^T & \mathbf{D}^T \\
\mathbf{C}^T & \mathbf{E}^T
\end{bmatrix}
$$

We transpose each matrix within the block matrix $\mathbf{H}$, and we also swap the upper and lower triangular matrices, i.e. in this case, $\mathbf{C}$ and $\mathbf{D}$. If we take a naive approach and set all matrices equal to $\mathbf{A}$, then we have for $\mathbf{H}$ and $\mathbf{H}^T$:

$$
\mathbf{H} =
\begin{bmatrix}
\mathbf{A} & \mathbf{A} \\
\mathbf{A} & \mathbf{A}
\end{bmatrix}\qquad
\mathbf{H}^T =
\begin{bmatrix}
\mathbf{A}^T & \mathbf{A}^T \\
\mathbf{A}^T & \mathbf{A}^T
\end{bmatrix}
$$

Since $\mathbf{A}$ is no longer symmetric, we cannot write $\mathbf{A}^T=\mathbf{A}$, and so, $\mathbf{H}^T\ne\mathbf{H}$. Let's try to remove the diagonal elements and see what happens:

$$
\mathbf{H} =
\begin{bmatrix}
0 & \mathbf{A} \\
\mathbf{A} & 0
\end{bmatrix}\qquad
\mathbf{H}^T =
\begin{bmatrix}
0 & \mathbf{A}^T \\
\mathbf{A}^T & 0
\end{bmatrix}
$$

Well, perhaps as expected, we still see that $\mathbf{H}^T\ne\mathbf{H}$. But we are getting closer, at least the diagonal entry looks good now. We have one final trick up our sleves. We know that the transpose of a transpose is just the orginal matrix itself, that is, $(\mathbf{A}^T)^T=\mathbf{A}$. So, let's write our lower triangular block matrix as the transpose of $\mathbf{A}$ and see what happens then:

$$
\mathbf{H} =
\begin{bmatrix}
0 & \mathbf{A} \\
\mathbf{A}^T & 0
\end{bmatrix}\qquad
\mathbf{H}^T =
\begin{bmatrix}
0 & (\mathbf{A}^T)^T \\
\mathbf{A}^T & 0
\end{bmatrix} =
\begin{bmatrix}
0 & \mathbf{A} \\
\mathbf{A}^T & 0
\end{bmatrix}
$$

Bingo! We have found a matrix $\mathbf{H}$ that is symmetric under transposition, i.e. we have $\mathbf{H}^T = \mathbf{H}$. So, we can apply our Conjugate Gradient algorithm now to $\mathbf{H}$, and it doesn't matter what $\mathbf{A}$ is, i.e. it can be symmetric or non-symmetric! Let's write our a full linear system for us to solve, i.e. $\mathbf{Hx}=\mathbf{b}$. We have:

$$
\begin{bmatrix}
0 & \mathbf{A} \\
\mathbf{A}^T & 0
\end{bmatrix}
\begin{bmatrix}
\tilde{\mathbf{x}}\\
\mathbf{x}
\end{bmatrix}=
\begin{bmatrix}
\mathbf{b} \\
\tilde{\mathbf{b}}
\end{bmatrix}
$$

Here, we see that we have two separate vectors for $\mathbf{x}$ and $\mathbf{b}$, where the second vector receives a tilde, i.e. $\tilde{\mathbf{x}}$ and $\tilde{\mathbf{b}}$. This equation now solves two separate problems, i.e. we can form two separate systems as:

$$
\mathbf{Ax}=\mathbf{b}\\
\mathbf{A}^T\tilde{\mathbf{x}}=\tilde{\mathbf{b}}
$$

We call this method Bi-Conjugate Gradient, with the prefix bi indicating that we have two systems to solve. The one with the tilde is also called the shadow system and our job now is to find a way to combine both as we compute a solution that will lead to convergence, despite $\mathbf{A}$ not being symmetrical. We can think of this as well like a predictor-corrector scheme; first, we compute the solution using the original system, i.e. $\mathbf{Ax}=\mathbf{b}$, which is then corrected by the shadow system $\mathbf{A}^T\tilde{\mathbf{x}}=\tilde{\mathbf{b}}$.

Let's start by writing out the residual vector at iteration $k$. For both systems we have:

$$
\mathbf{r}^k = \mathbf{b} - \mathbf{A}\mathbf{x}^k\\[1em]
\tilde{\mathbf{r}}^k =\tilde{\mathbf{b}} -\mathbf{A}^T\tilde{\mathbf{x}}^k
$$

Instead of enforcing orthogonality ($(\mathbf{r}^k)^T \mathbf{r}^{k+1}$) between subsequent iterations, we enforce bi-orthogonality between both the original and the shadow system, i.e. $\mathbf{r}$ and $\tilde{\mathbf{r}}$. We write this as:

$$
(\tilde{\mathbf{r}}^i)^T \mathbf{r}^j =0\quad \text{for} i \ne j
$$

Here, $i$ and $j$ can be though of as our iteration counter $k$, but we don't simply have $i=k$ and $j=k+1$. Instead, what this signifies is that all iterations $1,2,3,...,k-1$ are orthogonal to $k$. It is a much stronger orthogonality condition that simply requiring residuals at $k$ and $k+1$ to be orthogonal. To see, this, it is best to look at a figure that shows this condition. Let have a look at the following:

<!-- figure, width: 700px -->
![This figure shows how with each iteration k we build up the residuals of the shadow system that are all orthogonal to the residual of the original system at k+1](..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\BiCG_bi-orthogonality.png)

Here, I show the first three iterations of the Bi-Conjugate Gradient algorithm (with the first iteration being the initial guess of the residual vector $\tilde{\mathbf{r}}^0$). We can see that in the first iteration $k=0$ we have a simple line for our residual vector of the shadow system. We can't visualise the residual of the original system yet. Think of it as being a line in a different 1D system.

As we go to the next iteration $k=1$, we can see that $\mathbf{r}^1$ is orthogonal to $\tilde{\mathbf{r}}^0$. Both vectors form a plane in the 2D coordinate system. We can also see $\mathbf{r}^0$, i.e. the starting residual vector of the original system, and we see that it is not necessarily orthogonal to $\tilde{\mathbf{r}}^0$.

At the next iteration $k=2$, we can see that $\mathbf{r}^2$ is orthogonal to both $\tilde{\mathbf{r}}^0$ and $\tilde{\mathbf{r}}^1$. So, in this case, if we look at the bi-orthogonality equation where we had:

$$
(\tilde{\mathbf{r}}^i)^T \mathbf{r}^j =0\quad \text{for} i \ne j
$$

We can say that in the last iteration of $k=2$, we had $i=0,1$ and $j=2$. So, $\mathbf{r}^j$ is orthogonal to all previous shadow residuals $(\tilde{\mathbf{r}}^i)^T$. We also notice that with each iteration $k$ we take, we keep expanding the dimensionality of our system, and for that reason, I stopped at $k=2$, because I don't feel like attempting to draw $\mathbf{r}^3$ in 4D space. You'll have to trust me, though, that it exists, we just can't visualise it.

Similar to the bi-orthogonality, we also enforce this conditions for our search directions $\mathbf{d}$ as:

$$
(\tilde{\mathbf{d}}^i)^T\mathbf{A}\mathbf{d}^j = 0\quad\text{for } i\nej
$$

Let's pause for a moment, and take in what we have done, because if we understand this step, then the rest is just redoing the Conjugate Gradient algorithm derivation, which should be straightforward at this point. I like to look at analogies to help me understand things, so perhaps this one will help you grasp the concept of what the Bi-Conjugate Gradient method is doing.

Think of noise-cancelling earphones. We have a signal we want to hear (in our analogy, we have a system $\mathbf{Ax}=\mathbf{b}$ that we want to solve), but this signal contains noise from our surrounding environment. The chatter of people in a cafe, the comforting hum of an aircraft engine 10km above the pacific, or your annoying neighbour who decides that everyday is rock-o'clock.

In the absence of noise cancelling headphones, by the way, the correct solution to your annoying neighbour is primarily driven by speaker wattage, with your speakers directly attached to your ceiling so that your neighbour can fully enjoy the beat of [Helene Fischer](https://www.youtube.com/watch?v=haECT-SerHk). Not my prefered *beat*, but that was sort of the point. Let's just say I feel sorry for all of the other neighbours but I regret nothing!

Shortly after that back and forth, my neighbour had to report a broken water pipe and water was dripping into my flat (I wonder why ... If you live in the UK, it even rains inside your flat). I'm not sure if he had to pay for it, but let's just say he moved out shortly after that. I digress ...   

So, if the music that we want to hear with our headphones represents $\mathbf{Ax}=\mathbf{b}$, then the noise that we hear ontop of our music represents the residual $\mathbf{r} = \mathbf{b} - \mathbf{Ax}$. If we can cancel out all the noise, then our residual is zero, and we recover $\mathbf{Ax}=\mathbf{b}$ (or $\mathbf{b} - \mathbf{Ax}=0$).

If we want to cancel the noise, we need to create ani-noise (i.e. we record the noise in the surrounding and then flip the signal, hoping that this will cancel out all noise from our signal). So we have "music + noise + anti-noise = music" if you want, where "noise = - anti-noise". I know, highly scientific.

In our analogy, if we apply the Conjugate Gradient algorithm now to $\mathbf{Ax}=\mathbf{b}$, knowing that $\mathbf{A}$ is non-symmetric, we will make an error with each iteration, and that will lead to divergence. This error that we make in each iteration is the noise in our headphone analogy. We want to remove the noise (we want to correct the error of solving $\mathbf{Ax}=\mathbf{b}$), and so we introduce anti-noise (a shadow residual $\tilde{\mathbf{r}}^k$), which removes the noise (corrects the error introduced by the non-symmetric matrix).

The anti-noise has to have a specific structure, i.e. we can't jsut add any noise, it needs to be the exact opposite of the noise around us for noise-cancelling to work. In our analogy, we can't just take any residual, but it needs to be "the opposite" in some way. This opposite is expressed by the transposition of our coefficient matrix, i.e. $\mathbf{A}^T$.

It is related to $\mathbf{A}$ (the original system) so that if $\mathbf{A}$ changes, $\mathbf{A}^T$ will change as well and there is a certain connection between the two, just like the noise recorded by the microphone on our noise-cancelling headphones which is then flipped is related the the noise we hear. If the noise changes, the singla the microphone records will change in the same way. Thus, the noise and anti-noise are connected just like $\mathbf{A}$ and $\mathbf{A}^T$ are.

For any given time interval, our headphone will emit anti-noise ontop of what we are already hearing, and this is like taking one iteration in the Bi-Conjugate Gradient algorithm. This will then successfully cancel the noise, and this is like enforcing orthogonality in the Bi-Conjugate Gradient sense. I have summarised this comparison in the following table as well:

<!-- table, caption: "Comparison between Bi-Conjugate Gradient and noise-cancelling headphones" -->
| Headphones concept          | Bi-Conjugate Gradient meaning              |
| --------------------------- | ------------------------------------------ |
| Noise                       | Error / residual ( $\mathbf{r}^k$ )        |
| Anti-noise                  | Shadow residual ( $\tilde{\mathbf{r}}^k$ ) |
| Microphone measuring noise  | Action of ($A^T$)                          |
| Speaker emitting anti-noise | BiCG update step                           |
| Silence (cancellation)      | Orthogonality                              |

Hopefully this makes the Bi-Conjugate Gradient method somewhat clearer in what it does. As I said, if we are happy with this, then the rest of the algorithm will feel very natural. It is essentially the same as the Conjugate Gradient method, but we just take combinations from the original and shadow system to correct the errors introduced by the presence of the non-symmetrical matrix.

Let's start by reviewing the algorithm. Since the derivation is essentially the same as in the Conjugate Gradient method (and I must have grown a few more white hairs since I started to write this article 2 months ago), I'll focus on the differences instead of a full derivation. So, just to warn, this is the only section where I won't be deriving everything from scratch as it has already been done. 

We start with our line search update, i.e. finding a new location $\mathbf{x}$ on our surface $f$ that we want to minimise:

$$
\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k \mathbf{d}^k
$$

We found the residual to be:

$$
\mathbf{r}^{k+1} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k
$$

So far, so cosmopolitan. Now we make our first subtle change and also introduce the residual for the shadow system:

$$
\tilde{\mathbf{r}}^{k+1} = \tilde{\mathbf{r}}^k - \alpha^k \mathbf{A}^T\tilde{\mathbf{d}}^k
$$

The definition is exactly the same as the original system, except for the matrix $\mathbf{A}^T$. Remember, this is coming from our matrix $\mathbf{H}$, which was the 2-by-2 block matrix with only the matrix $\mathbf{A}$ and its transpose in it, so that $\mathbf{H}$ is symmetric.

The next step is to enforce orthogonality. We said that the residuals of the original system all have to be orthogonal to the shadow system's residuals at all previous iterations, though, for the algorithm itself, we only enforce orthogonality between two subsequent residuals or the original and shadow system:

$$
(\tilde{\mathbf{r}}^k)^T \mathbf{r}^{k+1} = 0
$$

We can insert the definition for $\mathbf{r}^{k+1}$ and get:

$$
(\tilde{\mathbf{r}}^k)^T(\mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k) = 0
$$

From this, we can find $\alpha^k$ as

$$
\alpha^k = \frac{(\tilde{\mathbf{r}}^k)^T\mathbf{r}^k}{(\tilde{\mathbf{r}}^k)^T\mathbf{A}\mathbf{d}^k}
$$

Let's focus on the denominator for a moment, i.e. $(\tilde{\mathbf{r}}^k)^T\mathbf{A}\mathbf{d}^k$. It is commonly written in a different form, so let's look at that. We know that the search direction can be updated (from the Conjugate Gradient algorithm) as:

$$
\mathbf{d}^{k+1} = \mathbf{r}^{k+1} + \beta^k \mathbf{d}^k
$$

Equally, the search direction of the shadow system can be updated as:

$$
\tilde{\mathbf{d}}^{k+1} = \tilde{\mathbf{r}}^{k+1} + \beta^k \tilde{\mathbf{d}}^k
$$

In both of these equations, we have the residual given at $k+1$, but for our denominator that we want to transform, we have the residual given at iteration $k$. We can write the direction update for the previous iteration as:

$$
\tilde{\mathbf{d}}^{k} = \tilde{\mathbf{r}}^{k} + \beta^k \tilde{\mathbf{d}}^{k-1}
$$

We just subtracted one from each iteration counter $k$, which resulted in the update to be based on the direction at $k-1$. Solving this for $\mathbf{r}^k$ we get:

$$
\tilde{\mathbf{r}}^{k} = \tilde{\mathbf{d}}^{k} - \beta^k \tilde{\mathbf{d}}^{k-1}
$$

Now we can insert that into our denominator and get:

$$
(\tilde{\mathbf{r}}^k)^T\mathbf{A}\mathbf{d}^k = (\tilde{\mathbf{d}}^{k} - \beta^k \tilde{\mathbf{d}}^{k-1})^T\mathbf{A}\mathbf{d}^k
$$

Let's expand the parenthesis:

$$
(\tilde{\mathbf{d}}^{k} - \beta^k \tilde{\mathbf{d}}^{k-1})^T\mathbf{A}\mathbf{d}^k = (\tilde{\mathbf{d}}^{k})^T\mathbf{A}\mathbf{d}^k - (\beta^k \tilde{\mathbf{d}}^{k-1})^T\mathbf{A}\mathbf{d}^k
$$

Now we have to remember that we enforce conjugancy. We see that we have $(\tilde{\mathbf{d}}^{k-1})^T\mathbf{A}\mathbf{d}^k$ in the second term, and two subsequent search directions. Conjugancy requires that two subsequent search directions multiplied by the matrix $\mathbf{A}$ is zero (this is the key difference between the Conjugate Gradient algorithm and the Steepest Descent method), thus, we can also write:

$$
(\tilde{\mathbf{d}}^{k} - \beta^k \tilde{\mathbf{d}}^{k-1})^T\mathbf{A}\mathbf{d}^k = (\tilde{\mathbf{d}}^{k})^T\mathbf{A}\mathbf{d}^k - \underbrace{(\beta^k \tilde{\mathbf{d}}^{k-1})^T\mathbf{A}\mathbf{d}^k}_{=0}\\[1em]
(\tilde{\mathbf{d}}^{k} - \beta^k \tilde{\mathbf{d}}^{k-1})^T\mathbf{A}\mathbf{d}^k = (\tilde{\mathbf{d}}^{k})^T\mathbf{A}\mathbf{d}^k
$$

Therefore, we have:

$$
(\tilde{\mathbf{r}}^k)^T\mathbf{A}\mathbf{d}^k = (\tilde{\mathbf{d}}^{k})^T\mathbf{A}\mathbf{d}^k
$$

With this, we can write $\alpha^k$ as:

$$
\alpha^k = \frac{(\tilde{\mathbf{r}}^k)^T\mathbf{r}^k}{(\tilde{\mathbf{r}}^k)^T\mathbf{A}\mathbf{d}^k} = \frac{(\tilde{\mathbf{r}}^k)^T\mathbf{r}^k}{(\tilde{\mathbf{d}}^k)^T\mathbf{A}\mathbf{d}^k}
$$

I have shown this step explicitly here as it is also usually skipped in the literature. This is a perfect point to pause for a second and mention that when you look up the Conjugate Gradient algorithm in the literature (and its derivatives, like the Bi-Conjugate Gradient algorithm), you will find subtle changes in the algorithms themselves. There are many equivalent ways of how we can write the algorithms, so stick with one source, ideally, and don't mix equations from separate sources. This *may* lead to an incorrect implementation.

For the search direction update, we already saw the definitions of $\mathbf{d}^{k+1}$ and $\tilde{\mathbf{d}}^{k+1}$. We want to use them now to construct $\beta$ (which will make subsequent search directions conjugate). We start again with orthogonality between the original and shadow system (a constraint coming from all the way back of the Steepest Descent method):

$$
(\tilde{\mathbf{r}}^{k+1})^T \mathbf{r}^k = 0
$$

We derived the equation for $\beta$ and that took a bit more space. If you do the same derivation now with $(\tilde{\mathbf{r}}^{k+1})^T$, instead of $(\mathbf{r}^{k+1})^T$, you'll arrive at:

$$
\beta^k = \frac{(\tilde{\mathbf{r}}^{k+1})^T \mathbf{r}^{k+1}}{(\tilde{\mathbf{r}}^k)^T\mathbf{r}^{k}}
$$

And this is the Bi-Conjugate Gradient algorithm! We need to make some choices about the initial vectors, and we have a few constraints. We have one constraint for the initial residual (at $k=0$), which is:

$$
(\tilde{\mathbf{r}}^{0})^T \mathbf{r}^0 \ne 0
$$

That is, residuals of the original and shadow system should not be the same orthogonal. We typically jsut set them equal at the first iteration,. i.e. $\tilde{\mathbf{r}}^{0} = \mathbf{r}^0$, which ensures they are co-linear (pointing in the same direction), which means they cannot be orthogonal.

The only other choice we make, similar to the Conjugate Gradient algorithm is that the first search direction is equivalent to the first residual for both systems, i.e.:

$$
\mathbf{d}^0 =\mathbf{r}^0,\quad \tilde{\mathbf{d}}^0 =\tilde{\mathbf{r}}^0
$$

The Bi-Conjugate Gradient algorithm has a tendency to oscillate around the minimum and thus it usually needs to be stabilised in order to be useful. As a result, we would't use the Bi-Conjugate Gradient algorithm method in practice. However, we do use the stabilised version of this algorithm, so let's have a look at this as well.

#### Summary of the Bi-Conjugate Gradient method

In summary, the Bi-Conjugate Gradient method can be written down as follows:

1. Compute initial residual: $\mathbf{r}^{0} = \mathbf{b} -\mathbf{Ax}^{0}$.
2. Set $\mathbf{d}^0 =\mathbf{r}^0,\quad \tilde{\mathbf{d}}^0 =\tilde{\mathbf{r}}^0 = \mathbf{r}^0$.
3. Loop from $k=0$ to $k=N$
4. Compute search direction length: $\alpha^k = \frac{(\tilde{\mathbf{r}}^k)^T\mathbf{r}^k}{(\tilde{\mathbf{d}}^k)^T\mathbf{A}\mathbf{d}^k}$
5. Compute new position on $f$: $\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k \mathbf{d}^k$.
6. Compute new residual: $\mathbf{r}^{k+1} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k$ and $\tilde{\mathbf{r}}^{k+1} = \tilde{\mathbf{r}}^k - \alpha^k \mathbf{A}^T\tilde{\mathbf{d}}^k$
7. Compute correction along $\mathbf{d}$: $\beta^k = \frac{(\tilde{\mathbf{r}}^{k+1})^T \mathbf{r}^{k+1}}{(\tilde{\mathbf{r}}^k)^T\mathbf{r}^{k}}$
8. Compute new (conjugate) search direction: $\mathbf{d}^{k+1} = \mathbf{r}^{k+1} + \beta^k \mathbf{d}^k$ and $\tilde{\mathbf{d}}^{k+1} = \tilde{\mathbf{r}}^{k+1} + \beta^k \tilde{\mathbf{d}}^k$
9. Compute solution residual: $res = \mathbf{x}^{k+1} - \mathbf{x}^{k}$ or $res = \mathbf{b} - \mathbf{Ax}^{k+1}$
10. If $||res||_2 \lt tolerance$ or $k=N$, break, otherwise, go back to step 3.
11. End

### The Bi-Conjugate Gradient Stabilised (BiCGStab) Gradient method

So, I said that the solution of the Bi-Conjugate Gradient method can be oscillatory. Does this sound familiar? When we went from the Steepest Descent method to the Conjugate Gradient algorithm, we said that the search directions in the Steepest Descent method introduced this oscillation, and the Conjugate Gradient algorithm got rid of it by introducing conjugate search directions, not orthogonal search directions.

So, where is the oscillation coming from in the Bi-Conjugate Gradient algorithm? Well, it isn't due to the search directions; these are still conjugate. But, since we use essentially the same (Conjugate Gradient) algorithm to compute the solution even though our matrix $\mathbf{A}$ may be non-symmetric, our residual may introduce oscillations. In the Bi-Conjugate Gradient algorithm, the residual was computed as:

$$
\mathbf{r}^{k+1} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k
$$

We introduced the shadow system for exactly that reason; get rid of any non-physical results introduced by this equation, i.e.:

$$
\tilde{\mathbf{r}}^{k+1} = \tilde{\mathbf{r}}^k - \alpha^k \mathbf{A}^T\tilde{\mathbf{d}}^k
$$

As it turns out, this will eventually converge, just like the Steepest Descent method will eventually converge, but it introduces oscillations which will slow down convergence. In fact, it can even lead to divergence in the worst case. So, we need to stabilse the residual calculation somehow. Let's see how we can achieve that.

There probably is one truth in CFD that we can use whenever we face oscillations: iteratively improve whatever is causing the oscillation. When we later solve the linear system of equations $\mathbf{Ax}=\mathbf{b}$ for non-linear equations, well, we iteratively apply $\mathbf{Ax}=\mathbf{b}$ until $\mathbf{x}^{n+1,k+1}\approx\mathbf{x}^{n+1,k}$, where $n$ is the time index and $k$ the iteration counter.

Laplacian smoothing will solve some form of laplacian equation to smooth a given field for a certain number of iterations $N$. For example, in mesh generation, we may smooth the grid point distributions to get a smoother grid to help us compute a smoother solution. In the residual smoothing technique we, well, smooth our residuals. This results in greater stability and explicit time stepping techniques, which are otherwise limited to a CFL number of 1, can be run up to a CFL number of 5 or so (personal experience).

In incompressible flows, we introduce smoothing (iterations) to couple the pressure to the velocity field. In fact, this is the main difference between the [SIMPLE](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-solve-incompressible-and-compressible-flows-in-cfd#aioseo-pressure-projection-based-simple) and the [PISO](https://cfd.university/learn/10-key-concepts-everyone-must-understand-in-cfd/how-to-solve-incompressible-and-compressible-flows-in-cfd#aioseo-how-to-go-from-the-simple-algorithm-to-the-piso-algorithm) algorithm, where just additional iterations are introduced to stabilise the pressure field.

Later in this article (yes, we are still not done, I'd suggest to put on the kettle and then throw all of your teabags in the kettle, this is not a one tea cup situation, you'll need the entire kettle), we will look at the multigrid approach, and in it, we will introduce smoothing to, you guessed it, stabilise the solution on each grid level.

We use the same idea here; instead of accepting the residual as it is, we are simply going to iterate over the residual, trying to minimise it in some regard. So, let's take the residual, which was given as:

$$
\mathbf{r}^{k+1} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k
\tag{eq:bicgstab-r-non-stabilised}
$$

And replace this with an intermediate residual:

$$
\mathbf{s}^{k} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k
\tag{eq:bicgstab-intermediate-residual}
$$

We now use this intermediate residual $\mathbf{s}$ and find a new residual $\mathbf{r}$ at $k+1$ as:

$$
\mathbf{r}^{k+1} = \mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k
\tag{eq:bicgstab-r-stabilised}
$$

Notice how Eq.(\ref{eq:bicgstab-r-stabilised}) retains the same form of Eq.(\ref{eq:bicgstab-r-non-stabilised}); we have just exchanged $\alpha^K$ here by $\omega^k$. So, you may say, that's stupid, we have created an intermediate residual $\mathbf{s}$, and the end result is that we have introduced an additional unknown $\omega$. But no, this isn't stupid; this is exactly what we want!

If we kept $\alpha^k$ in Eq.(\ref{eq:bicgstab-r-stabilised}), then all quantities on the right-hand side would be fixed. This means we would have no way of minimising $\mathbf{r}^{k+1}$. But, if we replace $\alpha^k$ by a new parameter $\omega^k$ of which we don't know the value, we can try to find a way to evaluate $\omega^k$ such that $\mathbf{r}^{k+1}$ is minimised. This is our goal.

Whenever we want to find the minimum of a function $f(x)$, we know that a possible minimum exists at $f'(x)=0$, i.e. the derivative is zero. It can be either a minimum, a maximum, or a saddle point. In our case, we want to find the minimum of our surface $f$, and we know that it only has a single minimum, with no maximum or saddle points, so we know that $f'=0$ will always be a minimum.

Now, looking at Eq.(\ref{eq:bicgstab-r-stabilised}) again, we see that the left-hand side is a vector, while the right-hand side only contains an unknown scalar $\omega^k$. We cannot minimise every entry in $\mathbf{r}^{k+1}$ if we only have a single scalar to tune. Thus, we try to minimise some appropriate norm of $\mathbf{r}^{k+1}$ instead.

Taking the norm of a vector, e.g. the $L_2$-norm $||\mathbf{r}^{k+1}||_2$, has the property that a vector is reduced to a scalar value. If we then minimise this scalar value, we know that this will overall minimise the underlying vector, i.e. in this case $\mathbf{r}^{k+1}$. Since we want to minimise the norm of $\mathbf{r}^{k+1}$, we have to take the norm on both the left-hand side and right-hand side in Eq.(\ref{eq:bicgstab-r-stabilised}). This results in:

$$
||\mathbf{r}^{k+1}||_2 = ||\mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k||_2
\tag{eq:bicgstab-r-stabilised-norm}
$$

So, we are left with the task of minimising the right-hand side, i.e.:

$$
\text{min}||\mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k||_2
$$

In general, for a vector $\mathbf{v}$, the $L_2$-norm is defined as:

$$
||\mathbf{v}||_2 = \sqrt{v_1^2 + v_2^2 + v_3^2 + ... + v_N^2}
$$

Using Eq.(\ref{eq:bicgstab-r-stabilised-norm}), we can write this as:

$$
||\mathbf{r}^{k+1}||_2 = \sqrt{(s_1^k - \omega^k As_1^k)^2 + (s_2^k - \omega^k As_2^k)^2 + (s_3^k - \omega^k As_3^k)^2 + ... + (s_N^k - \omega^k As_N^k)^2}
$$

Here, $As_1^k$, $As_2^k$, $As_3^k$, etc. are the products of the entire vector $\mathbf{s}^k$ with the first row, second row, third row, etc. of $\mathbf{A}$. Now we have the parameter $\omega^k$ within the square-root, so we need to square the expression to be able to solve for $\omega^k$. Squaring results in:

$$
||\mathbf{r}^{k+1}||_2^2 = ||\mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k||_2^2
\tag{eq:bicgstab-r-stabilised-norm}
$$

This can be written as:

$$
||\mathbf{r}^{k+1}||_2 = (s_1^k - \omega^k As_1^k)^2 + (s_2^k - \omega^k As_2^k)^2 + (s_3^k - \omega^k As_3^k)^2 + ... + (s_N^k - \omega^k As_N^k)^2
$$

We could also write this as:

$$
||\mathbf{r}^{k+1}||_2 = (s_1^k - \omega^k As_1^k)\cdot(s_1^k - \omega^k As_1^k) + (s_2^k - \omega^k As_2^k)\cdot(s_2^k - \omega^k As_2^k) + (s_3^k - \omega^k As_3^k)\cdot(s_3^k - \omega^k As_3^k) + ... + (s_N^k - \omega^k As_N^k)\cdot(s_N^k - \omega^k As_N^k)
$$

This is nothing else than the scalar (dot) product, and so, we can write this with vector quantities as:

$$
||\mathbf{r}^{k+1}||_2^2 = (\mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k)\cdot(\mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k)
$$

We can also write this as a multiplication of a row with a colum vector:

$$
||\mathbf{r}^{k+1}||_2^2 = (\mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k)^T(\mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k)
$$

This is what we looked at already previously, i.e. we saw that we can write $\mathbf{a}\cdot\mathbf{b}=\mathbf{a}^T\mathbf{b}$ when we derived the Conjugate Gradient algorithm. Using this equation, we can multiply terms to get get rid of the parenthesis:

$$
||\mathbf{r}^{k+1}||_2^2 = (\mathbf{s}^k)^T\mathbf{s}^k - \omega^k(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k - \omega^k(\mathbf{A}\mathbf{s}^k)^T\mathbf{s}^k + (\omega^k)^2(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k
\tag{eq:bicgstab-omega-expanded}
$$

The rules of transposition state that for two vectors $\mathbf{u}$ and $\mathbf{v}$, the following rule holds:

$$
\mathbf{u}^T\mathbf{v} = \mathbf{v}^T\mathbf{u}
$$

Let's verify this with the following example:

$$
\mathbf{u}=
\begin{bmatrix}
a\\b\\c
\end{bmatrix},\qquad
\mathbf{v}=
\begin{bmatrix}
d\\e\\f
\end{bmatrix}
$$

We have:

$$
\mathbf{u}^T\mathbf{v}=
\begin{bmatrix}
a,b,c
\end{bmatrix}
\begin{bmatrix}
d\\e\\f
\end{bmatrix}=
ad+be+cf
$$

We also have:

$$
\mathbf{v}^T\mathbf{u}=
\begin{bmatrix}
d,e,f
\end{bmatrix}
\begin{bmatrix}
a\\b\\c
\end{bmatrix}=
da+eb+fc=ad+be+cf
$$

This shows that the equality $\mathbf{u}^T\mathbf{v} = \mathbf{v}^T\mathbf{u}$ does indeed hold for this example (and, in general, for any other real vectors). Using this equality, and setting $\mathbf{u}=(\mathbf{s}^k)$ and $\mathbf{v}=\mathbf{As}^k$, we can see that in Eq.(\ref{eq:bicgstab-omega-expanded}), we have these two terms:

$$
-\omega^k(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k - \omega^k(\mathbf{A}\mathbf{s}^k)^T\mathbf{s}^k
$$

We can also write this as:

$$
-\omega^k \mathbf{u}^T\mathbf{v} - \omega^k \mathbf{v}^T\mathbf{u}
$$

Using $\mathbf{u}^T\mathbf{v} = \mathbf{v}^T\mathbf{u}$, we can rewrite the second term as:

$$
-\omega^k \mathbf{u}^T\mathbf{v} - \omega^k \mathbf{u}^T\mathbf{v}
$$

This can be simplified to:

$$
-2\omega^k \mathbf{u}^T\mathbf{v}
$$

Substituting $\mathbf{u}$ and $\mathbf{v}$ back into our equations, we get:

$$
-\omega^k(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k - \omega^k(\mathbf{A}\mathbf{s}^k)^T\mathbf{s}^k = -2\omega^k(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k
$$

We can use this to simplify Eq.(\ref{eq:bicgstab-omega-expanded}) to:

$$
||\mathbf{r}^{k+1}||_2^2 = (\mathbf{s}^k)^T\mathbf{s}^k - \omega^k(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k - \omega^k(\mathbf{A}\mathbf{s}^k)^T\mathbf{s}^k + (\omega^k)^2(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k\\[1em]
||\mathbf{r}^{k+1}||_2^2 = (\mathbf{s}^k)^T\mathbf{s}^k - 2\omega^k(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k + (\omega^k)^2(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k
\tag{eq:bicgstab-omega-expanded}
$$

We differentiate this equation now with respect to $\omega$:

$$
\frac{\partial ||\mathbf{r}^{k+1}||_2^2}{\partial \omega} = - 2(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k + 2(\omega^k)(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k
$$

If we want to minimise this function now, we have to set the derivative (left-hand side) to zero and then solve for $\omega$. We get:

$$
0 = - 2(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k + 2(\omega^k)(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k\\[1em]
2(\omega^k)(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k = 2(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k\\[1em]
\omega^k(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k = (\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k\\[1em]
\omega^k = \frac{(\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k}{(\mathbf{A}\mathbf{s}^k)^T\mathbf{A}\mathbf{s}^k}
$$

And this is it. We have found $\omega^k$ for iteration $k$, so we can evaluate it now so that we are guaranteed to get a new residual $\mathbf{r}^{k+1}$ which is minimised. This will avoid spurious oscillations once we are close to the minimum. Our job is to now include this new parameter $\omega^k$ in our derivation, so that we avoid said oscillations.

In the Bi-Conjugate Gradient algorithm, we updated the new location on the surface $f$ as:

$$
\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k \mathbf{d}^k
$$

This corresponded to a residual of:

$$
\mathbf{r}^{k+1} = \mathbf{b}-\mathbf{Ax}^{k+1}
$$

However, in the Bi-Conjugate Gradient Stabilised version, we changed the definition for the residual $\mathbf{r}$, by introducing the intermediate residual $\mathbf{s}$, which resulted in the additional parameter $\omega^k$. So, we need to use this new definition for the residual and then reverse engineer, if you want, the update equation for $\mathbf{x}^{k+1}$. Let's do that.

To get started, let's introduce the residual $\mathbf{r}^{k+1/2}$. This is an intermediate residual, and you can think of it as the residual we would obtain from the Bi-Conjugate Gradient algorithm, i.e. with no correction applied to it to minimise $\mathbf{r}$. Going from $\mathbf{r}^{k+1/2}$ to $\mathbf{r}^{k+1}$ is what we do in the stabilised version of the Bi-Conjugate Gradient algorithm, i.e. our residual correction.

We can write this intermediate residual as:

$$
\mathbf{r}^{k+1/2} = \mathbf{b}-\mathbf{Ax}^{k+1/2}
$$

Equally, we can find the location before we apply any correction to the residual as:

$$
\mathbf{x}^{k+1/2} = \mathbf{x}^k + \alpha^k \mathbf{d}^k
$$

Inserting $\mathbf{x}^{k+1/2}$ into the equation for $\mathbf{r}^{k+1/2}$, we get:

$$
\mathbf{r}^{k+1/2} = \mathbf{b}-\mathbf{Ax}^{k+1/2}\\[1em]
\mathbf{r}^{k+1/2} = \mathbf{b}-\mathbf{A}(\mathbf{x}^k + \alpha^k \mathbf{d}^k)\\[1em]
\mathbf{r}^{k+1/2} = \underbrace{\mathbf{b}-\mathbf{A}\mathbf{x}^k}_{\mathbf{r}^k} - \alpha^k \mathbf{A}\mathbf{d}^k\\[1em]
\mathbf{r}^{k+1/2} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k
$$

By comparing this with Eq.(\ref{eq:bicgstab-intermediate-residual}), which was given as:

$$
\mathbf{s}^{k} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k
$$

We see that we have found $\mathbf{r}^{k+1/2}=\mathbf{s}^k$. Now we return to the difinition for the residual definition of the Bi-Conjugate Gradient Stabilised algorithm, which was given in Eq.(\ref{eq:bicgstab-r-stabilised}) as:

$$
\mathbf{r}^{k+1} = \mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k
$$

This is one equation for $\mathbf{r}^{k+1}$, but we also know the general definition for the residual i.e. $\mathbf{r}^{k+1} = \mathbf{b}-\mathbf{Ax}^{k+1}$. We can set both of these equal and get:

$$
\mathbf{b}-\mathbf{Ax}^{k+1} = \mathbf{s}^k - \omega^k \mathbf{A}\mathbf{s}^k
$$

In this equation, we see that we have $\mathbf{x}^{k+1}$ given on the left-hand side, and we want to isolate that. Once we have done that, we know what our update equation is for the new position $\mathbf{x}^{k+1}$ on our surface $f$. First, we replace $\mathbf{s}^k$ with $\mathbf{r}^{k+1/2}$. This gives us:

$$
\mathbf{b}-\mathbf{Ax}^{k+1} = \mathbf{r}^{k+1/2} - \omega^k \mathbf{A}\mathbf{s}^k
$$

We know that we have two expressions for $\mathbf{r}^{k+1/2}$. We compute it using either of thw following two expressions:

$$
\mathbf{r}^{k+1/2} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k\\[1em]
\mathbf{r}^{k+1/2} = \mathbf{b} - \mathbf{Ax}^{k+1/2}
$$

Let's see what happens when we insert the first expression. We get:

$$
\mathbf{b}-\mathbf{Ax}^{k+1} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k - \omega^k \mathbf{A}\mathbf{s}^k
$$

We can subtract $\mathbf{b}$ and multiply by $-1$ and get:

$$
\mathbf{Ax}^{k+1} = -\mathbf{r}^k + \alpha^k \mathbf{A}\mathbf{d}^k + \omega^k \mathbf{A}\mathbf{s}^k + \mathbf{b}
$$

Since we want to isolate $\mathbf{x}^{k+1}$ on the left-hand side, we need to multiply by the inverse of $\mathbf{A}$, i.e. $\mathbf{A}^{-1}$ and get:

$$
\mathbf{A}^{-1}\mathbf{Ax}^{k+1} = -\mathbf{A}^{-1}\mathbf{r}^k + \alpha^k \mathbf{A}^{-1}\mathbf{A}\mathbf{d}^k + \omega^k \mathbf{A}^{-1}\mathbf{A}\mathbf{s}^k + \mathbf{A}^{-1}\mathbf{b}\\[1em]
\mathbf{I}\mathbf{x}^{k+1} = -\mathbf{A}^{-1}\mathbf{r}^k + \alpha^k \mathbf{I}\mathbf{d}^k + \omega^k \mathbf{I}\mathbf{s}^k + \mathbf{A}^{-1}\mathbf{b}
$$

Since a vector multiplied by the identity matrix is the vector itself, we can remove that from the calculation. Don't believe me? Well, I think you and me need to work on your trust issue, but just to show that I am not using black magic here, take the following example:

$$
\begin{bmatrix}
1 & 0 & 0\\
0 & 1 & 0\\
0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
a \\ b \\ c
\end{bmatrix} =
\begin{bmatrix}
1\cdot a + 0\cdot b + 0\cdot c\\
0\cdot a + 1\cdot b + 0\cdot c\\
0\cdot a + 0\cdot b + 1\cdot c
\end{bmatrix} =
\begin{bmatrix}
a + 0 + 0\\
0 + b + 0\\
0 + 0 + c
\end{bmatrix} =
\begin{bmatrix}
a\\
b\\
c
\end{bmatrix}
$$

Happy? Then let's remove the identity matrix from our equation. We get:

$$
\mathbf{x}^{k+1} = -\mathbf{A}^{-1}\mathbf{r}^k + \alpha^k \mathbf{d}^k + \omega^k \mathbf{s}^k + \mathbf{A}^{-1}\mathbf{b}
$$

We can now insert $\mathbf{r}^k=\mathbf{b}-\mathbf{Ax}^k$ and get:

$$
\mathbf{x}^{k+1} = -\mathbf{A}^{-1}(\mathbf{b}-\mathbf{Ax}^k) + \alpha^k \mathbf{d}^k + \omega^k \mathbf{s}^k + \mathbf{A}^{-1}\mathbf{b}\\[1em]
\mathbf{x}^{k+1} = -\mathbf{A}^{-1}\mathbf{b} + \mathbf{A}^{-1}\mathbf{Ax}^k + \alpha^k \mathbf{d}^k + \omega^k \mathbf{s}^k + \mathbf{A}^{-1}\mathbf{b}\\[1em]
\mathbf{x}^{k+1} = -\mathbf{A}^{-1}\mathbf{b} + \mathbf{I}\mathbf{x}^k + \alpha^k \mathbf{d}^k + \omega^k \mathbf{s}^k + \mathbf{A}^{-1}\mathbf{b}\\[1em]
\mathbf{x}^{k+1} = -\mathbf{A}^{-1}\mathbf{b} + \mathbf{x}^k + \alpha^k \mathbf{d}^k + \omega^k \mathbf{s}^k + \mathbf{A}^{-1}\mathbf{b}\\[1em]
$$

We have both $-\mathbf{A}^{-1}\mathbf{b}$ and $+\mathbf{A}^{-1}\mathbf{b}$ in this equation, so we can cancel it. This results in:

$$
\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k \mathbf{d}^k + \omega^k \mathbf{s}^k
$$

We have found an update equation, and, if we compare it against the update equation for $\mathbf{x}^{k+1}$ from the Bi-Conjugate Gradient algorithm, we had:

$$
\mathbf{x}^{k+1} = \mathbf{x}^k + \alpha^k \mathbf{d}^k
$$

We see that the stabilised version introduces the additional term $\omega^k \mathbf{s}^k$. Without it, we get oscillatory residuals, with it, we smooth those residuals. I said that we have two equations for $\mathbf{r}^{k+1/2}$ that we coudl have used, i.e. $\mathbf{r}^{k+1/2} = \mathbf{r}^k - \alpha^k \mathbf{A}\mathbf{d}^k$ (this is what we used) and $\mathbf{r}^{k+1/2} = \mathbf{b} - \mathbf{Ax}^{k+1/2}$.

If we insert the second definition, that is, $\mathbf{r}^{k+1/2} = \mathbf{b} - \mathbf{Ax}^{k+1/2}$, we will obtain exactly the same result.

Ok, so we have a new position on the surface $f$, next we need to figure out what our search direction is. Let's remind ourselves what $\beta$ was, as shown in the following figure which we saw previously:

<!-- wp:image {"width":"400px","sizeSlug":"large","align":"center"} -->
<figure class="wp-block-image aligncenter size-large is-resized"><img src="..\..\articles\07_10-key-concepts-everyone-must-understand-in-cfd\assets\05_5_matrix\conjugate_gradient_new_direction.png" alt="A sketch showing how the new search direction d at k+1 points along what was the error e at k+1, which shows that the conjugate gradient method points towards the minimum and thus can achieve faster convergence than the Steepest Descent method." class="wp-image-5550" style="width:400px"/></figure>
<!-- /wp:image -->

$\beta^k$ is the correction we apply so that subsequent search directions are conjugate to each other (or A-orthogonal). Since the residual is, at least initially, the same as the search direction $\mathbf{d}$, we see that changes to the residuals also changes the search directions, and so we need to incorporate our correction factor $\omega^k$ into our update equation for $\beta^k$ as well.


















This will ensure that we are not starting to oscillate on $f$ close to the minimum. We also modify the computation of $\beta$ from:

$$
\beta^k = \frac{(\tilde{\mathbf{r}}^0)^T \mathbf{r}^{k+1}}{(\tilde{\mathbf{r}}^0)^T \mathbf{r}^k}
$$

in the Bi-Conjugate Gradient method to:

$$
\beta^k = \frac{(\tilde{\mathbf{r}}^0)^T \mathbf{r}^{k+1}}{(\tilde{\mathbf{r}}^0)^T \mathbf{r}^k} \cdot \frac{\alpha^k}{\omega^k}
$$

in the Bi-Conjugate Gradient Stabilised version, where this value is now scaled as a fraction of both $\alpha^k$ and $\omega^k$. Finally, a new search direction is found as:

$$
\mathbf{d}^{k+1} = \mathbf{r}^{k+1} + \beta^k \left( \mathbf{d}^k - \omega^k A \mathbf{d}^k \right)
$$

Here, we subtract $\omega^k A \mathbf{d}^k$ from $\mathbf{d}^k$ in the stabilised version, again, to reduce oscillations.


### The Generalised Minimal Residual (GMRES) method

## Preconditioning

### The Jacobi preconditioner

### The Incomplete Cholesky (IC) precoinditioner

### The Incomplete Lower-Upper (ILU) preconditioner

## The multigrid approach

### Geometric multigrid approaches

### Algebraic multigrid approaches

### Multigrids as preconditioners

## Extension of linear system of equations solver to non-linear systems

### The Picard-iteration

### The Newton-Raphson approach

## Summary