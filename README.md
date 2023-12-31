
# vesicle-quantification

Quantification of vesicle size and deformation in TEM images.

## Approach

These scripts use cryosparc 4 (Structura Biotechnology, https://structura.bio/) as a user interface to define structures of membrane vesicles in negative stain or cryo TEM images.


## Dependencies

The data processing scripts require python 3.5 or higher, preferably installed in a dedicated virtual environment.  The scripts should be placed into a dedicated working directory. Package dependencies are:
 - numpy
 - json
 - math
 - statistics
 - pillow
 - scipy
 - matplotlib

Both `cs_to_spline_interpolation_cs4.py` and `halir_flusser_ellipse.py` need to be stored in your working directory.

In addition to the python environment, you will also need a system with [proc2d](https://blake.bcm.edu/emanwiki/EMAN1/Programs/Proc2D) installed.

## Preprocess micrographs and select particles in cryosparc

 - Load your images (movies of still images) into cryosparc.
 - Perform patch motion correction and patch ctf estimation on the imported images.  Motion correction needs to be performed even if the images are single-frame.
   - Take note of your motion correction job identifier.
 - Load your images into a manual picking job.
 - Set particle size to something small (~10 Angstrom) and choose evenly spaced particles around the outside of your vesicles.

![image](https://github.com/tribell4310/vesicle-quantification/assets/67428134/b2592a9c-9bff-4a47-8d25-f7b08d32d8c9)


 - When annotating multiple vesicles on a single micrograph, annotate one on one side, then the other.
   - Downstream scripts will look for jumps in how far sequential picks are apart to discriminate beween vesicles, so try to maximize how far apart you can make picks between sequential vesicles!
 - When finished, click "Done Picking" to terminate the picking job.

## Retrieve defined coordinate info from cryosparc

 - Navigate in your filesystem browser to your cryosparc project directory followed by the folder for your manual picking job.
 - Within this folder there will be a file called `extracted_particles.cs`.  Copy this file to your local working directory.

## Convert preprocessed images

 - Navigate in your filesystem browser to your cryosparc project directory followed by the folder for your patch motion correction job.
 - Within this folder, enter the `motioncorrected` subdirectory.
 - Within your local working directory, create two new folders called `png_out2` and `splines`.
 - Copy all files from `motioncorrected` that end with `_patch_aligned_doseweighted.mrc` to `png_out2`.
 - In your local `png_out2` folder, run the `convert2.sh` script:

`bash convert2.sh`

 - This will convert the mrc files to png images.
 - Delete the input mrc files, leaving only the png images in your `png_out2` folder.

## Analyze particles with python

- In an appropriate python environment, run the vesicle fitting and quantification script.

`python cs_to_spline_interpolation_cs4.py motionCorrJobID cryosparcParticleFile`

- The input requires your motioncorr job ID - for example, if your motioncorr job ID was J6 and your particle file was `extracted_particles.cs`, then the command would be:

`python cs_to_spline_interpolation_cs4.py J6 extracted_particles.cs`

 - Images of vesicle fits to your micrographs can be found in the `splines` subdirectory.  Quantification of the vesicles can be found in the `Vesicle_data` subdirectory.
   - Vesicle area
   - Vesicle perimiter
   - Vesicle effective circular diameter
   - Vesicle circularity (1 is a perfect circle)
   - Vesicle stratification order (for multilamellar samples, how many other vesicles contain this vesicle)

![016843626466719142702_prom_wt_bucket_all_020_patch_aligned_doseweighted_SpLiNeS](https://github.com/tribell4310/vesicle-quantification/assets/67428134/f5096912-84b7-4a73-9d41-3abf4817b337)


![017404409993085251713_prom_wt_bucket_all_044_patch_aligned_doseweighted_SpLiNeS](https://github.com/tribell4310/vesicle-quantification/assets/67428134/1e4cc8b6-f761-47a3-a255-4741bdae8865)

![image](https://github.com/tribell4310/vesicle-quantification/assets/67428134/91d0d331-4d17-4544-a970-8bc3e4a50f63)



## Questions

If you run into issues with these scripts, please feel free to open an Issue in GitHub, and I'll get back to you as quickly as I can.  Thanks for your interest!


> Written with [StackEdit](https://stackedit.io/).
