# Cookies Documentation

For this package, we recommend using cookiecutter (along with Flask or Dash which is built with Flask)
to set up your directory if you intend to create a web app/interface for your project.

`Cookies` makes it very easy to do this.

Learn more about the [cookiecutter](https://github.com/audreyr/cookiecutter) package.

## Examples

The [Manager module](https://github.com/datasnakes/OrthoEvolution/tree/master/OrthoEvol/Manager)
uses the _CookBook_ and _Oven_ classes as a primary means of functioning.

### A Simple Implementation

```python
from OrthoEvol.Cookies import Oven

Kitchen = Oven(repo="repo", user="username", project="project-name",
               output_dir="path/to/project")
Pantry = Kitchen.Ingredients

# To create a user directory
Kitchen.bake_the_user()
```
