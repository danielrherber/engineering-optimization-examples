# README (engineering-optimization-examples)

The repository contains several examples for helping teach optimization concepts.
These are utilized in a course on Engineering Optimization ([ENGR 510](https://www.online.colostate.edu/courses/ENGR/ENGR510.dot)) taught by [Dr. Daniel R. Herber](https://github.com/danielrherber) at Colorado State University.
All examples are available in Matlab, and many have Python equivalents.

## Summary

| Topic        | Folder          | # of Examples |
|--------------|-----------------|---------------|
| Introduction | [/1-introduction](1-introduction) | 4        |
| Linear Optimization | [/2-linear-optimization](2-linear-optimization) | 4        |
| Unconstrained Optimization | [/3-unconstrained-optimization](3-unconstrained-optimization) | 5         |
| Constrained Optimization | [/4-constrained-optimization](4-constrained-optimization) | 4         |
| Derivative-free Optimization | [/5-derivative-free-optimization](5-derivative-free-optimization) | 4          |

## Running the Matlab Examples

1. Install an appropriate Matlab version and toolboxes to your machine
	- Validated version is ``R2022a`` but likely compatible with recent older versions
	- *[CSU-only installation instructions]* [https://www.engr.colostate.edu/ets/matlab/](https://www.engr.colostate.edu/ets/matlab/)
	- Required toolboxes:
		- Symbolic Math Toolbox
		- Optimization Toolbox
		- Global Optimization Toolbox
1. Open an example of interest (often in a ``/*/matlab/`` directory) and run the example
	- Examples might have some combination of generated command window text and figures 
1. *[Optional]* Run [test_matlab_examples.m](test_matlab_examples.m) to verify that all examples work
	- You will need to make sure all project files are in your path
	- There are some required user inputs as well (simple ``Enter`` commands will suffice)

## Running the Python Examples

1. Install an appropriate python version to your machine
	- Validated version is ``3.9.x`` but likely compatible with recent older versions
	- [https://www.python.org/downloads/](https://www.python.org/downloads/)
1. Install the following packages:
	- ``numpy``
	- ``scipy``
	- ``matplotlib``
	- ``sympy``
	- ``networkx``
1. Navigate to an example of interest (often in a ``/*/python/`` directory) and run the example
	- Examples might have some combination of generated text output and figures 