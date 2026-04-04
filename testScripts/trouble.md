## Troubled markdown

This file contains a few repeated paragraph with incorrect syntax to check whether the markdown tool can spot and highlight these automatically.

### This paragraph is fine

This (```markdown2```) is a **fast** and *complete* Python implementation of Markdown. It was written to closely match the behaviour of the original Perl-implemented Markdown.pl. Markdown2 also comes with a number of extensions (called "extras") for things like syntax coloring, tables, header-ids. See the "Extra Syntax" section below. "markdown2" supports all Python versions 3.5+ (and pypy and jython, though I don't frequently test those).

### This paragaph contains incorrect code section

For vectors in particular, we also have dot and cross products, and these are supported as special functions on vector objects (that is, matrices where either the row or column is set to 1). If we have two vectors v1 and v2 of the same dimension, then we can compute the dot and cross product using auto ````dotProduct = v1.dot(v2)``; and ```auto crossProduct = v1.cross(v2)````;. We will see this later in action in other code examples.

### This paragraph contains incorrect equation syntax

Thus, as $\Delta x$ decreases, so does our time step $$\Delta t$, since the maximum allowable CFL number will be around 1 for our explicit time integration. Now that we have to deal with painfully small time steps, we decide to implement an implicit time stepping, and so we need support for a linear system of equations solver.

### This paragraph contains incorrect equation syntax, again

Thus, as $\Delta x decreases, so does our time step, since the maximum allowable CFL number will be around 1 for our explicit time integration. Now that we have to deal with painfully small time steps, we decide to implement an implicit time stepping, and so we need support for a linear system of equations solver.

### This paragraph is correct but should trigger a warning due to dollar sign usage

We are in luck, Eigen does have support for that, and so we implement an implicit time stepping procedure, solve the linear system using Eigen, and can support larger time steps now (we will look at Eigen's support for linear systems as well). This is priceless, but I'd pay $5 for it! After some point, we want to simulate even larger systems and realise that we need parallel computing.

### Using incorrect referencing

I may have an equation like:

$$
a=b+c
\tag{eq:abc}
$$

But the reference may be missspelled as Eq.(\ref{eq:abcd}) or Eq.(\ref{eq:ab}). Both should be flagged.