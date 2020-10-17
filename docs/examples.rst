========
Examples
========

For examples, please ref `example` folder.

Note: For CASTEP, the user should copy `RunCASTEP.bat` or `RunCASTEP.sh` in their Materils Studio installation folder.

For VASP, the user should provide the `POTCAR` file

After the user provide the above files, just run `runElastic3rd.py` file by `python runElastic3rd.py`

If the user want to save the log into file, please run `python runElastic3rd.py >> Result.txt`

If the user want to run Elastic3rd in queue system, the following is an example. Please note that it assume the python and vasp is exists in the `PATH`


.. code::

    #PBS -l nodes=1:ppn=24
    #PBS -l walltime=48:00:00
    #PBS -l pmem=8gb
    #PBS -A open
    #PBS -j oe

    cd $PBS_O_WORKDIR
     
    python runElastic3rd.py