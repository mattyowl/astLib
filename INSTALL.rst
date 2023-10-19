astLib requires:

* Python (tested on versions 3.6+)
* Astropy - http://www.astropy.org (tested on version 3.2.1)
* Numpy - http://numpy.scipy.org (tested on version 1.18.1)
* SciPy - http://scipy.org (tested on version 1.3.1)
* Matplotlib - http://matplotlib.sourceforge.net (tested on version 3.1.1)

Optional:
   
* Python Imaging Library - http://www.pythonware.com/products/pil (tested on version 1.1.7)

Other versions of the software listed above are likely to work.

You can install astLib via pip:

.. code-block::

   pip install astLib --user


You may also install using the standard ``setup.py`` script, e.g., as root:

.. code-block::

   sudo python setup.py install


Alternatively, 

.. code-block::

   python setup.py install --user


will install ``astLib`` under ``$HOME/.local`` (on Ubuntu), and in some other default location on Mac.

You can also use the ``--prefix`` option, e.g.,

.. code-block::

   python setup.py install --prefix=$HOME/local


and then add, e.g., ``$HOME/local/lib/python3.6/site-packages`` to 
$PYTHONPATH (adjust the path according to your Python version number).

.. code-block::

   export PYTHONPATH=$HOME/local/lib/python3.6/site-packages:$PYTHONPATH


Installation on recent versions of macOS (may no longer be relevant)
====================================================================

Some users have reported that the standard method for installing ``astLib`` does not work on recent versions
of macOS (e.g., Big Sur), due to the default compiler flags. The current workaround for this is to install
using:
  
.. code-block::

   CFLAGS="-Wno-error=implicit-function-declaration" python setup.py install
   

Thanks to Michael Cowley and Stefano Covino for helping to resolve this issue.
