# mp_flux
This is a basic particle detection algorithm that was created as a side project to my PhD learning julia.

It is probably not the most efficient algorithm, but it is sufficient for my needs.

# how to use

Use an IDE of your choice, I really liked to use atom for this.

Source the basic_function.jl, source the mp.flux.jl

I uploaded three test files.

# Description of the algorithm

The first step is to determine a color intensity threshold (CIT), which converts the image into a bit image. Thus, either a pixel has a color intensity above the threshold and is set to one or it is below the threshold and it is set to zero.

The remaining pixels are checked for connections. Only horizontal and vertical connections are recognized. An
option to also recognize diagonal connections it not implemented.

The algorithm can return several outputs.

The standard outputs are two tables:

- particle_hist_data_[CIT].csv This table contains the size (number of pixels a particle contains) and the time (timestamp of the image the data was extracted from). For now there is also a type variable, which says "pre", which is a relict from older versions. It does not have a real meaning.
- setup_data_info_[CIT].csv This table contains information about the camera setup. 

The file name is saved with the used CIT value, so if you used multiple CITs on the same data set you can match the tables to the CIT.

Very usefull for me was:

- filter_image : Creates an grayscale image of the particles that were detected by the algorithm.

# image metadata

To use the image metadata functions you need to install exiftool.

https://exiftool.org/

# Julia version

I ran this code on Julia Version 1.5.3.


***This project is licensed under the terms of the MIT license.***
