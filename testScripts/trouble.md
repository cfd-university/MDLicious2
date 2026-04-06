## Troubled markdown

This file contains a few repeated paragraph with incorrect syntax to check whether the markdown tool can spot and highlight these automatically.

### This paragraph is fine

This (```markdown2```) is a **fast** and *complete* Python implementation of Markdown. It was written to closely match the behaviour of the original Perl-implemented Markdown.pl. Markdown2 also comes with a number of extensions (called "extras") for things like syntax coloring, tables, header-ids. See the "Extra Syntax" section below. "markdown2" supports all Python versions 3.5+ (and pypy and jython, though I don't frequently test those).

### This paragaph contains incorrect code section

For vectors in particular, we also have dot and cross products, and these are supported as special functions on vector objects (that is, matrices where either the row or column is set to 1). If we have two vectors v1 and v2 of the same dimension, then we can compute the dot and cross product using auto ````dotProduct = v1.dot(v2)``; and ```auto crossProduct = v1.cross(v2)````;. We will see this later in action in other code examples.

But we can also have code environments:

<!-- code, caption: "Wrong tag should be flagged, should start with tag{code:", \tag{cod:listing-1} -->
```python
print('hello world')
```

This should not be triggered, neither code sections without any identifier

```
>>> Hello, World
```

But this certainly shoudl fail:

``
wrong backticks detected
```

### This paragraph contains incorrect equation syntax

Thus, as $\Delta x$ decreases, so does our time step $$\Delta t$, since the maximum allowable CFL number will be around 1 for our explicit time integration. Now that we have to deal with painfully small time steps, we decide to implement an implicit time stepping, and so we need support for a linear system of equations solver.

### This paragraph contains incorrect equation syntax, again

Thus, as $\Delta x decreases, so does our time step, since the maximum allowable CFL number will be around 1 for our explicit time integration. Now that we have to deal with painfully small time steps, we decide to implement an implicit time stepping, and so we need support for a linear system of equations solver.

### Yet more equations:

Dangerous, ```markdown2``` may mess with LaTeX equations. It is fine in equation environments:

$$
\int_{0}^1 a\mathrm{d}x \approx \sum_{k=1}^N a_i\Delta x
$$

But dangerous inline, i.e. $\int_{0}^1 a\mathrm{d}x \approx \sum_{k=1}^N a_i\Delta x$

### This paragraph is correct but should trigger a warning due to dollar sign usage

We are in luck, Eigen does have support for that, and so we implement an implicit time stepping procedure, solve the linear system using Eigen, and can support larger time steps now (we will look at Eigen's support for linear systems as well). This is priceless, but I'd pay $5 for it! After some point, we want to simulate even larger systems and realise that we need parallel computing.

### Using incorrect referencing

I may have an equation like:

$$
a=b+c
\tag{eq:abc}
$$

Each equation should have their own tag, even if none is provided. The equation above shoudl be Eq.(\ref{eq:abc}), below no tag is provided:

$$
1+1=2
$$

One with a tag

$$
2+2=4
\tag{eq:obvs}
$$

And now one without again

$$
4+4=8
$$

All should get a tag, regardless.

Reference may be missspelled as Eq.(\ref{eq:abcd}) or Eq.(\ref{eq:ab}). Both should be flagged.
