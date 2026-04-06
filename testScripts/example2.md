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

Let's discretise that with a simple first-order Euler method in time, and a second-order central scheme in space. We will denote the unknown (next) time level by $n+1$, while all variables at the current time level $n$ are known. Therefore, we obtain the following equations: See Eq.(\ref{eq:unsteady-heat-diffusion})

- **Explicit** time integration:

$$
\frac{T_i^{n+1}-T_i^{n}}{\Delta t} = \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}
\tag{eq:explicit-heat-diffusion}
$$

- **Implicit** time integration: See Eq.(\ref{eq:explicit-heat-diffusion})

$$
\frac{T_i^{n+1}-T_i^{n}}{\Delta t} = \Gamma\frac{T^{n+1}_{i+1} - 2T^{n+1}_{i} + T^{n+1}_{i-1}}{(\Delta x)^2}
$$

See Eq.(\ref{eq:explicit-heat-diffusion})

We see that the only difference is at which time level we evaluate the temperatures on the right-hand side, i.e. either at $n$ or at $n+1$. The next step is to collect all unknowns on the left-hand side of the equation, all knowns on the right-hand side of the equation, and then to write the equation in coefficient form. This is done in the following way for both schemes:

- **Explicit** time integration:

$$
\frac{T_i^{n+1}}{\Delta t} = \frac{T_i^{n}}{\Delta t} + \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}\\[1em]
\left[\frac{1}{\Delta t}\right]T_i^{n+1} = \frac{T_i^{n}}{\Delta t} + \Gamma\frac{T^{n}_{i+1} - 2T^{n}_{i} + T^{n}_{i-1}}{(\Delta x)^2}
\tag{eq:explicit-system-discretised}
$$

See Eq.(\ref{eq:explicit-system-discretised})

- **Implicit** time integration:

$$
\frac{T_i^{n+1}}{\Delta t} - \Gamma\frac{T^{n+1}_{i+1} - 2T^{n+1}_{i} + T^{n+1}_{i-1}}{(\Delta x)^2} = \frac{T_i^{n}}{\Delta t}\\[1em]
\left[\frac{1}{\Delta t}\right]T_i^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i+1}^{n+1} + \left[\frac{2\Gamma}{(\Delta x)^2}\right]T_{i}^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i-1}^{n+1} = \frac{T_i^{n}}{\Delta t}\\[1em]
\left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i+1}^{n+1} + \left[\frac{1}{\Delta t} + \frac{2\Gamma}{(\Delta x)^2}\right]T_{i}^{n+1} + \left[\frac{-\Gamma}{(\Delta x)^2}\right]T_{i-1}^{n+1} = \frac{T_i^{n}}{\Delta t}
\tag{eq:implicit-system-discretised}
$$

See Eq.(\ref{eq:implicit-system-discretised})

A common notation here is to introduce coefficients $a_P$, $a_E$, and $a_W$ for the coefficients in-front of $T^{n+1}_{i}$, $T^{n+1}_{i+1}$, and $T^{n+1}_{i-1}$, respectively, with $E$ and $W$ standing for east and west. Furthermore, all of the known quantities on the right-hand side can be put together into a single value, which we will denote as $b$. Doing so gives us the following system of equations: