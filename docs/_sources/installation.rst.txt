Installation
============

There are two ways to install the ``battery_sizing_model`` package.

1. From PyPI through ``pip``.
2. From GitHub so you can modify the source code.

PyPI
----

To install from PyPI through ``pip`` simply run:

.. code-block:: console

    pip install battery_sizing_model

GitHub
------

If you want to download the whole repo to make changes to the source code and have
access to the data processing notebooks and scripts, clone or fork the 
`repo <https://github.com/seanbuchanan-eng/battery_sizing_model>`_
then run ``pip install -r dev-requirements.txt`` in a virtual environment running python 
3.10+.

Once all of the dependencies are installed, open ``data_processing/save_batch_to_disk.py``,
change the file path in line 98

.. code-block:: python

    raws_prepath = 'C:/Users/seanb/OneDrive/Documents/PRIMED/export/batdat/MTC/raws/'

to the path to data on your machine, and run the file. Running the file will process the data
for each cell separately and store it to disk before starting the processing of the next cell.
This ensures that the program doesn't run out of RAM. Upon completion, there should be a new 
``data_objects`` directory in your ``battery_sizing_model`` directory.

Finally, open ``data_processing/deg_average_model_data.py`` and run it. Once complete, there
should be a new ``aggregate_params`` directory in your ``battery_sizing_model`` directory.

After completing these steps all of the files in the repository should work.