# EPG
Extended Phase Graphs

This set of python 3.6 code is derived from Brian Hargreaves EPG code<sup>1</sup> written in Matlab. The main routines have been written in C++ to speed up the code and uses SWIG<sup>2</sup> to wrap the code and produce python functions. 

## Installation

The python files are created using the command 

```python setup.py build```

and if no errors have occured, then to install

```python setup.py install```


---

1. http://web.stanford.edu/~bah/software/epg
2. http://swig.org
