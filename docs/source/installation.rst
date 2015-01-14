.. _installation:

******************************
Prerequisites and installation
******************************

RNAchrom source code, installation scripts, and source code of this website are available at

	`RNAchrom github page`_

.. _nucChIP github page: https://github.com/mrivas/RNAchrom

Prerequisites
=============

To use RNAchrom, you need at least version 2.5 of `Python <http://www.python.org/>`_.


Installation on Linux
=====================

RNAchrom installation package will take care of all dependencies. You just need to download the *source* files, unzip them, and go into the unzipped directory.

.. code-block:: bash

   $ wget https://github.com/mrivas/RNAchrom/archive/master.zip RNAchrom.zip
   $ unzip RNAchrom.zip
   $ cd RNAchrom

To install RNAchrom locally ( only for the current user in case you don't have super-user permission) type

.. code-block:: bash

   $ python setup.py install --user

To make RNAchrom available to all users, use instead

.. code-block:: bash

   $ python setup.py build
   $ sudo python setup.py install

To test the installation

.. code-block:: bash

   $ python
   >> import RNAchrom

if no errors are reported, your are all set.
