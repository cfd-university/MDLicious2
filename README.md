![Version number](https://img.shields.io/badge/Version-0.24.1-red.svg)

## Introduction

MDLicious2 is a ```markdown``` to ```html``` converter based on [python-markdown2](https://github.com/trentm/python-markdown2) with additional conversion of additional (unsupported) Markdown syntax as used in [articles](https://github.com/cfd-university/articles) on [cfd.university](https://cfd.university).

The new syntax that is added allows to easily integrate:
- Figures, captions, and figure numbering
- Tables, captions, and table numbering
- Code listings, captions, and listing numbering
- Equations and its numbering
- Embedding of YouTube videos with privacy-enhanced (non-tracking) ```iframe```s

All, except the YouTube ```iframe```s, provide numbering support, and this is extended with the ability to reference each element with a corresponding ```\ref{...}``` label anywhere in the text. This is very much inspired by LaTeX and provides the same functionality. 

## Usage

To use this converter, you need to provide a json configuration file of the form:

```json
{
    "inputFile": "/path/to/markdown/file.md",
    "outputDirectory": "/output/directory/",
    "replace": {
        "find this string in text": "replace it with this text",
        "some more matches to replace": "with yet another string"
    }
}
```

Only the ```inputFile``` and ```outputDirectory``` are strictly speaking required. The generated HTML file will take the same filename as the markdown file, with only the extension changed from ```md``` to ```html```. Other inputs are optional. If this ```JSON``` file is stored as config.json in the current directory, then we can run it as:

```bash
python3 /path/to/MDLicious2.py /path/to/config.json
```

## Prerequisites

This converter requires Python3 to run, as well as the following packages:

- Markdown 2
- Pygments
- Beautiful Soup 4

To install these, create a virtual environment first:

```bash
python3 -m venv .venv
```

Then, source the virtual environment

```bash
# UNIX
. .venv/bin/active

# Windows
.venv/Scripts/Activate.ps1
```

You can now install these packages into the virtual environment with

```bash
pip install -r requirements.txt
```

In addition, LaTeX equations are directly converted into HTML through the Node.js module katex. Therefore, Node.js and npm needs to be installed, as well as katex itself. Download [Node.js](https://nodejs.org/en/download) and [npm](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm) and then install katex with the following command:

```bash
npm install katex
```

Verify that katex is indeed installed by running:

```bash
npx katex --version
```

If this prints a version number, you are probably good to go.

## New Markdown components

As alluded to in the introduction, new components have been added to Markdown. The syntax of each new component is described in the following sections.

### Figures

Markdow supports the use of figures using the following syntax:

```markdown
![Alt text of image](/source/of/image)
```

This allows to specify the alt text and the source to the image, but it has no mechanism for us to set the size of the image, a caption, or a tag to use later for referencing back to this image. To remedy this, an HTML comment is introduced above each figure like so:

```markdown
<!-- figure, width: 600px, caption: "This is the figure caption, which can be different from the alt text", \tag{fig:figure-1} -->
![Alt text of image](https://placehold.co/600x400)
```

The initial ```<!-- figure``` is what is used to identify this as a figure, so you will need to match that exactly. We are now able to specify the width of the image, provide a caption, and a tag. Using the tag, we can reference back to the figure. For example, from anywhere in a text, we may refer back to the figure using the ```\ref{fig:figure-1}``` label. Internally, this will be converted to a number.

> [!TIP]  
> The caption and ```\tag``` argument are optional. If you provide a caption, you will get something like "Figure 1: caption goes here". If you don't provide a caption, you simply get "Figure 1". If you do not provide a ```\tag```, one will be generated for you. You only need to provide it if you want to reference back to the figure later.

> [!CAUTION]
> When you do provide a tag, you must start it with ```\tag{fig:...}```. Change ... to the tag you want to use. This is important because any tag starting with ```fig:``` is interpreted as a figure.

### Tables

Tables follow a very similar syntax. We can define standard tables in ```markdown``` as:

```markdown
|        | Group 1 |     | Group 2      |     |
|--------|---------|-----|--------------|-----|
|        | A       | B   | A            | B   |
| Test 1 | 1.0     | 2.7 | 1.4          | 1.2 |
| Test 2 | 1.1     | 1.4 | 1.7          | 1.3 |
```

We can extend that with a comment just as we saw with the figure. An example is shown below:

```markdown
<!-- table, caption: "The table's caption", \tag{tab:table-1} -->
|        | Group 1 |     | Group 2      |     |
|--------|---------|-----|--------------|-----|
|        | A       | B   | A            | B   |
| Test 1 | 1.0     | 2.7 | 1.4          | 1.2 |
| Test 2 | 1.1     | 1.4 | 1.7          | 1.3 |
```

The initial ```<!-- table``` is what is used to identify this as a table, so you will need to match that exactly. We can also use the tag again in the text to reference back to this particular table with ```\ref{tab:table-1}```

> [!TIP]  
> The caption and ```\tag``` argument are optional. If you provide a caption, you will get something like "Table 1: caption goes here". If you don't provide a caption, you simply get "Table 1". If you do not provide a ```\tag```, one will be generated for you. You only need to provide it if you want to reference back to the table later.

> [!CAUTION]
> When you do provide a tag, you must start it with ```\tag{tab:...}```. Change ... to the tag you want to use. This is important because any tag starting with ```tab:``` is interpreted as a table.

### Code listings

Markdown supports code as well natively, but just like tables and figures, it is extended to allow for custom captions and tagging, but in a slightly different way as before. We can augment any code snippet with a code comment like so:

~~~markdown
<!-- code, caption: "caption to go under code listing", \tag{code:listing-1} -->
```python
print('hello cfd.university')
```
~~~

This will print the code listing with its caption underneath. However, since I sometimes want to use code sections to simply print results, bash comments, or the like, I don't want to have an automatic caption like "Listing 1", "Listing 2", "Listing 3", and so on.

The following logic is used:

- If a caption is present as an HTML comment, a caption will be extracted and printed.
- If a caption is not present, but a ```\tag``` is present, a simple "Listing 1", "Listing 2", "Listing 3", ... will be printed.
- If no HTML comment is present, nothing will be printed.

Examples:

~~~markdown
<!-- code, caption: "Hello world in Python", \tag{code:listing-1} -->
```python
print('hello cfd.university')
```
~~~

This will result in:

```python
print('hello cfd.university')
```
Listing 1: Hello world in Python

We can remove the caption and have:

~~~markdown
<!-- code, \tag{code:listing-1} -->
```python
print('hello cfd.university')
```
~~~

This will result in:

```python
print('hello cfd.university')
```
Listing 1

Or, we just remove the entire HTML meta data comment, as seen in the following:

~~~markdown
```python
print('hello cfd.university')
```
~~~

This will result in:

```python
print('hello cfd.university')
```

That's it.

> [!CAUTION]
> When you do provide a tag, you must start it with ```\tag{code:...}```. Change ... to the tag you want to use. This is important because any tag starting with ```code:``` is interpreted as a code listing.

### Equations

There are two separate types of equations we can use: inline and multiline. Inline equations appear within paragraphs. We insert equations using a leading and trailing dollar sign. So, for example, an inline equation such as $\mathbf{F}=m\mathbf{a}$ may be written as ```$\mathbf{F}=m\mathbf{a}$```.

You can use standard LaTeX equation syntax in your equations. For multiline equations, we start on a new line with two dollar signs, followed by the equations, and we finish with two dollar signs on a new line. For example, the following multiline equation:

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

can be written as:

```latex
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
```

Multiline equation can also contain a ```\tag``` if we want to reference equations later. ```\tag```s are written on a new line, just before the end (trailing two dollar signs), for example:

```latex
$$
\mathbf{F}=m\mathbf{a}
\tag{eq:newtons-second-law}
$$
```

Internally, this ```\tag``` will be replaced by a number. If you do not provide a ```\tag```, one will be automatically generated for you, so that equations are consecutively numbered. We can then reference these equations as ```Eq.(\ref{eq:newtons-second-law})```, for example. Multiline equation do not necessarily have to have more than one line. They will always be shown on a new line in their own environment.

> [!CAUTION]
> When you do provide a tag, you must start it with ```\tag{eq:...}```. Change ... to the tag you want to use. This is important because any tag starting with ```eq:``` is interpreted as an equation.

### Youtube video embedding

This is a feature not natively available in Markdown. In order to embbed a youtube video, simply provide an HTML comment as ```<!-- youtube: URL>```. You will need to change the ```URL``` here to the youtube video's URL. It can contain any number of variables in the URL, those will be ignored. For example, the following two YouTube video embbedings are equivalent:

```markdown
Display YouTube video:

<!-- youtube: https://www.youtube.com/watch?v=3gdWqgBrERk&t=831s -->

The same video but now without any arguments/variables:

<!-- youtube: https://www.youtube.com/watch?v=3gdWqgBrERk -->
```

Appropriate classes are automatically added so that the video scales to the correct size on the HTML page. Videos are embedded using a privacy-enhanced mode, that is, YouTube is not allowed to set any cookies on [cfd.university](https://cfd.university).

## MDLicious2 - What happend to version 1?

We shall not speak of version 1, ever, again, ...
