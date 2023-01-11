==================================
Waston Fitting & Tracking Examples
==================================

Below a number of directly executable commands for Watson fitting and tracking are listed. The relative paths are set so that the execution is done from the scripts folder. See Watson_README.rst for details on the scripts.

watson-fitting
~~~~~~~~~~~~~~

.. code-block:: console
    
    python watson-fitting --i ../../test_data/792766/fodf/mtdeconv/ankele/ --init given --initfile ../../test_data/792766/lowrank/fodf_peaks.nrrd -ob ../../outputfolder/watson_backup.npz -o ../../outputfolder/watson_tracking_data.nrrd

watson-backup-to-data
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: console
    
    python watson-backup-to-data --i ../../test_data/792766/watson/o4/watson_backup.npz -o ../../outputfolder/watson_tracking_data.nrrd -vvi ../../outputfolder/watson_vvi_cone_data


watson-tracking
~~~~~~~~~~~~~~~

Tracking with pre-computed Watson fit:

.. code-block:: console
    
    python prob-tracking --i ../../test_data/792766/fodf/mtdeconv/ankele/ --infile ../../test_data/792766/watson/o4/watson_tracking_data.nrrd --seedpoints ../../test_data/792766/seeding/CG_right/seeds_5000.pts -o ../../outputfolder/CG_right_5000.tck --prob Watson --interpolation FACT --wmmin 0.4 --rank 3

Tracking with fODF Interpolation:

.. code-block:: console
    
    python prob-tracking --i ../../test_data/792766/fodf/sh_fodf/o4/ --infile ../../test_data/792766/watson/o4/watson_tracking_data.nrrd --wmvolume ../../test_data/792766/fodf/mtdeconv/ankele/wmvolume.nrrd --seedpoints ../../test_data/792766/seeding/CG_right/seeds_5000.pts -o ../../outputfolder/CG_right_5000_intp.tck --prob Watson --interpolation TrilinearFODFWatson --wmmin 0.4 --rank 3 --dist 0 --var 6 --exp 3
