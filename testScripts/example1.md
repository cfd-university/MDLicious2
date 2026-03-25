Provide a table of content first:

[toc]

# Heading 1

## Heading 2

### Heading 3

#### Heading 4

##### Heading 5

###### Heading 6

This (```markdown2```) is a **fast** and *complete* Python implementation of Markdown. It was written to closely match the behaviour of the original Perl-implemented Markdown.pl. Markdown2 also comes with a number of extensions (called "extras") for things like syntax coloring, tables, header-ids. See the "Extra Syntax" section below. "markdown2" supports all Python versions 3.5+ (and pypy and jython, though I don't frequently test those).

## Testing links

I want to be able to embed links like [google.com](https://www.google.com) directly into my text.

## Equations

I want equations to render natively using katex, e.g. $\mathbf{F}=m\mathbf{a}$, my favourite example, should render in text, not using katex script tags. I may have more than one equation as $\mathbf{Ax}=\mathbf{b}$. Equations on new lines should work as well:

$$
\mathbf{F}=m\mathbf{a}
\tag{eq:newton}
$$

But also multiline equations:

$$
\begin{bmatrix}
1 \\
2 \\
3
\end{bmatrix}\cdot
\begin{bmatrix}
11 \\
22 \\
33
\end{bmatrix}
\tag{eq:matrix-equation}
$$

Of course, Eq.(\ref{eq:newton}) and Eq.(\ref{eq:matrix-equation}) can also be referenced now!

## Images

I want to be able to embed images using the native markdown syntax, on a new line, like so:

<!-- figure, width: 600px, caption: "This image is reproduced from [Teschner et al.](https://www.sciencedirect.com), I can even have more than one link, see: [google.com](https://www.google.com), fantastic! Though, I may also have equations to write $1 + 2 = 3$.", \tag{fig:figure_1} -->
![This is the alt text for an image, notice the comment above for extra information](https://placehold.co/600x400)

<!-- figure, width: 300px, caption: "Another figure", \tag{fig:figure-2} -->
![This is the alt text for an image, notice the comment above for extra information](https://placehold.co/300x100)

And now, I can refer to the Figure as Figure \ref{fig:figure_1} and Figure \ref{fig:figure-2}.

## Block quotes

I also want block quotes to display correctly. For example

> Clever stuff someone once said
>
> [Clever dick](https://www.wikipedia.org)

This should show as indented

## Lists

Ah, yes, lists, someone else should take care of this. I need single list, enumerated and itemised, but also nested lists, so let's test:

- item 1
- item 2
- item 3

This should be no problem. Neither should be:

1. Item 1
2. Item 2
3. Item 3

How about nested lists?

- Heading 1
- Heading 2
  - Heading 2.1
  - Heading 2.2
    - Heading 2.2.1
  - Heading 2.3
- Heading 3

## Code

Another fancy issue. Code in text should resolve natively, like ```code``` should be formated differently. But, here is something for pygments to take care of:

<!-- code, caption: "A code section with a subtitle", \tag{code:python-hello-world} -->
```python
import numpy as np

def main():
    print('Hello numpy')
```

I should be able now to refer back to the code as Listing \ref{code:python-hello-world}. Also, I want to have support for C++. What happens if I don't provide any caption?

```c++
#include <iostream>

int main(int *argc, char[] *argv) {

  // important string
  std::cout << "Hello C++" << std::endl;
  return 0;
}
```

Now another code section with a caption:

<!-- code, caption: "Just another code section ...", \tag{code:python-hello-world-2} -->
```python
import numpy as np

def main():
    print('Hello numpy')
```

This should be Listing 2, we can veryfy that by refering back to its tag (Listing \ref{code:python-hello-world-2}). Is it the same? Brilliant, that means that the second code sectionw as ignored!

<!-- code, \tag{code:python-hello-world-3} -->
```python
import numpy as np

def main():
    print('Hello numpy')
```

A code section without a ```caption``` but with a ```tag``` should still yield a section beneath the code with "Listing 1", "Listing 2", "Listing 3", and so on. I should be able to refer back to it as well, for example, see Listing \ref{code:python-hello-world-3}.


## Youtube

I also want to embed youtube video directly into my code, without having to think about it, and assign the correct classes or styles.

<!-- youtube: https://www.youtube.com/watch?v=3gdWqgBrERk&t=831s -->

This should convert to:

## Tables

What would markdown conversion be without support for tables? Here is one:

<!-- table, caption: "This caption should appear above the table from now on. And yes, captions here should also support equations $\mathbf{Ax}=\mathbf{b}$!", \tag{tab:table-1} -->
|        | Group 1 |     | Group 2      |     |
|--------|---------|-----|--------------|-----|
|        | $A$     | B   | $\mathbf{A}$ | B   |
| Test 1 | 1       | 2   | 14           | 423 |
| Test 2 | text    | Tom | more te4xt   |     |

Now let's print the same table again, but this time, I won't use any captions:

<!-- table -->
|        | Group 1 |     | Group 2      |     |
|--------|---------|-----|--------------|-----|
|        | $A$     | B   | $\mathbf{A}$ | B   |
| Test 1 | 1       | 2   | 14           | 423 |
| Test 2 | text    | Tom | more te4xt   |     |

The above should print a table! Can I also refer back to Table \ref{tab:table-1}? Can I refer back several time? Table \ref{tab:table-1}, Table \ref{tab:table-1}, Table \ref{tab:table-1}!