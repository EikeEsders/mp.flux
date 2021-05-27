# mp.flux
This is a basic particle detection algorithm that was created as a side project to my PhD learning julia.

It is probably not the most efficient algorithm, but it worked so far for me ;)

# What can it be used for ?

I used it to separate particles from a background. One first step is to determine a color intensity threshold (CIT),
which converts the image into a bit image. Thus, either a pixel has a color intensity above the threshold and
is set to one or it is below the threshold and it is set to zero.

The remaining pixels are checked for connections. Only horizontal and vertical connections are recognized. An
option to also recognize diagonal connections it not implemented.

The algorithm can return several outputs.

The standard outputs are two tables:

- particle_hist_data_[CIT].csv This table contains the size (number of pixels a particle contains) and the time (timestamp of the image the data was extracted from). For now there is also a type variable, which says "pre", which is a relict from older versions. It does not have a real meaning.
- setup_data_info_[CIT].csv This table contains information about the camera setup.

The file name is saved with the used CIT value, so if you used multiple CITs on the same data set you can match the tables to the CIT.

# Julia version

I ran this code on Julia Version 1.5.3.


***This project is licensed under the terms of the MIT license.***
