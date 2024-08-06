# mypackage/__init__.py

# This will allow users to import gene_density and graphs directly from the package
from .gene_density import *
from .graphs import *

# Optionally, you can define a list of modules to be exported when `import *` is used on the package
__all__ = ['gene_density', 'graphs']