====================================
RMG - Reaction Mechanism Generator
====================================

.. image:: https://github.com/GreenGroup/RMG-Java/raw/master/web/source/_static/rmg_logo.png
   :align: right
      
RMG is an automatic chemical reaction mechanism generator that constructs
kinetic models composed of elementary chemical reaction steps using a general
understanding of how molecules react.

The easiest way to get started on Windows is to download the complete installer package
from over at http://rmg.sourceforge.net/, and read the documentation_ that's also hosted there.
If you want to use the cutting-edge development version that is hosted here at GitHub,
then you'll have to compile it yourself.

.. _documentation: http://rmg.sourceforge.net/documentation/

Compiling from source
-----------------------------

Linux
~~~~~

Mac
~~~~~

You will need the Java compiler, Apache ant, and GNU Make, but these are all provided by MacOS.
If you're missing these `developer tools <http://developer.apple.com/technologies/tools/>`_ you can download them from Apple.
You will also need a Fortran compiler. 
We recommend gfortran_, but g95 can be used if you add the option `F90=g95` to the `make` command. 
On MacOS X 10.6 (Snow Leopard) the gfortran compiler flags can be optimized with the make option `MACOS=true`.
One nice way to get gfortran is to install homebrew_ then then type `brew install gfortran`.

.. _gfortran: http://r.research.att.com/tools/
.. _homebrew: http://mxcl.github.com/homebrew/

The quick-start instructions in full:

#. git clone git://github.com/GreenGroup/RMG-Java.git
#. cd RMG-Java
#. make MACOS=true
#. make test

For more, refer to the documentation_.


Windows
~~~~~~~~
