from setuptools import setup
from Cython.Build import cythonize

setup(
    name="gravity_noise",
    ext_modules=cythonize(
        "gravity_noise.pyx",
        compiler_directives={"language_level": "3"}
    ),
    zip_safe=False,
)
