
mp_flux("2021_03_29_75_125_B/"
        , 0.3
        , filter_image = true
        , area_density = true
        , filter_color = "red"
        , from = 1
        , to = 1
        , N_f = 128
        , write_results = false)







## threshold evaluation

original_crop_png("2021_03_08_215_B/threshold_evaluation/2021_03_08_215_B01321.tif"
        , 5896, 5896 + 830 , 2522 , 2522 + 605)

threshold_evaluation("2021_03_08_215_B/threshold_evaluation/"
        , 15    # value counted visually
        , 0.04    # lowest threshold
        , 0.98    # highest threshold
        , 0.02   # interval raised every iteration
        , "red"
        , 5896, 5896 + 830 , 2522 , 2522 + 605)

#=

NECESSARY

dir = data directory
intensity_filter::Float64 = 0-1; threshold for selecting pixels in the green light spectrum
filter_color::string = "green"
vertical_steps::Int64 = 10; number of vertical subsets the image is devided into; 5 is enough
timestamp_mode::String = "metadata", "postifx" ; either use the metadata of images or a postfix to write a time information into the result table


OPTIONAL USEFULL

intensity_filter_mode::String = "absolute", "relative"; absolute relates the threshold to the absolute maximum of a green value 255 in RGB; relative relates the threshold to the max green value found in the image (changes from image to image);
dont_write_results::Bool = false, true; if set to false there will be now results written into a csv file, this might be interesting if you just want to compute filter or marked images
use_kernel:Bool = true; uses a 3x3 mean kernel to smooth the image before the intensity filter; useful if there is noise in your particles, but they are distinct from there suroundings; I used it at ISO 51200


SUBSET ANALYSIS

from::Int64 = 1 ;  set the first image to be analysed
to::Int64 = 0 ;  set the last image to be analysed

crop::Bool = false, ; if set true, you can analyse only a fraction of the image; you  have to set the boundaries of the crop; only rectangles are possible
        x1::Int64 = 1,
        x2::Int64 = 100,
        y1::Int64 = 1,
        y2::Int64 = 100


MORE DETAILED PARTICLE ANALYSIS

circularity_analysis::Bool = false ; adds the circularity of each particle to the result table
image_format::String = "JPG", "tif"
aspect_ratio_mode::Bool = false ;  adds the aspect ratio of each particle to the result table
position_calculation::Bool = false ;  adds the y and x position of each particle to the result table

DETAILED IMAGE ANALYSIS

intensity_filter_image::Bool = false ; output an png image of all pixels that were above the threshold
timestamp_mode::String = "metadata", "postifx" ; either use the metadata of images or a postfix to write a time information into the result table
mark_mode::Bool = false ; output an png image of all registered particles
mark_size::Int64 = 40  ;  set the mark size, length in pixels

=#
