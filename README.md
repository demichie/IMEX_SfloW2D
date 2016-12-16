## Welcome to IMEX_SfloW2D


IMEX_SfloW is a FORTRAN90 code designed to solve the shallow water equations, written as a hyperbolic system of partial differential equations with relaxation and source terms. The model is discretized in time with an explicit-implicit Runge-Kutta method where the hyperbolic part is solved explicetely and the other terms (relaxation and surce) are treated implicitely. The finite volume solver for the hyperbolic part of the system is based on a semidiscrete central scheme and it is not tied on the specific eigenstructure of the model. The implicit part is solved with a Newton-Raphson method where the elements of the Jacobian of the nonlinear system are evaluated numerically with a complex step derivative technique.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/demichie/IMEX_SfloW2D/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
