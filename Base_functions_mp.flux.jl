using Images, Statistics, DelimitedFiles, StatsBase, ImageFiltering, OffsetArrays

function center_of_mass_address(particles, N_f, pic_height)

    N_p = length.(particles[:])  # convert to length of every array (N_p)

    n = findall(x -> x >= N_f, N_p) # define filter and get all arrays with more values than N_f

    N_fi = particles[n] #  index all arrays with length \beq than N_f

    #### reduction to 1 pixel size ####

    lin_address = zeros(Int64,length(N_fi),3)

    for i = 1:length(N_fi) # calculate center of mass

        y = mod1.(N_fi[i], pic_height)

        x = cld.(N_fi[i], pic_height)

        y = round(sum(y) / length(y), digits = 0)

        x = round(sum(x) / length(x), digits = 0)

        lin_address[i,1] = x

        lin_address[i,2] = y

        lin_address[i,3] = 1
    end

    return lin_address

end

function re_image_center_of_mass(pic_width, pic_height, particles, center,dir, file)

    img = zeros(pic_width * pic_height) # creates an "empty" bit-array

    for s = 1:length(particles)
        color = 0.3

        for i in particles[s]
            img[i] = color[1]

        end

    end

    for i in center
        img[i] = 0.8

    end

    img_post = reshape(img, pic_height, :)

    save(string(dir, "/",file,"center_of_mass.png"),img_post
    )

end

function threshold_evaluation(dir, by_eye, from, to, interval, filter_color,N_f, x1 ,x2 ,y1 ,y2)

        by_eye = by_eye

        values = collect(from:interval:to)

        threshold_limit = Int64[]

        for i in values

                try

                mp_flux(dir, i, image_format = "tif",filter_image = true, crop = true, x1 = x1, x2 = x2, y1 = y1, y2 = y2, use_kernel = false,N_f = N_f, filter_color = filter_color)

                catch info

                if isa(info, BoundsError)

                        threshold_limit = vcat(threshold_limit,i)

                elseif isa(info, ArgumentError)

                        threshold_limit = vcat(threshold_limit,i)

                end

                end

        end

        files = readdir(dir, join = true)

        filter!(s -> occursin(r"hist",s), files)

        result = hcat("#","threshold")

        for i in files

                particles = size(readdlm(i))[1]-1

                threshold = "0." * split(SubString(split(i,"_")[end]),".")[2] # not sure what happenend here...

                entry = hcat(particles, threshold)

                result = vcat(result, entry)

        end

        a = result[2:end,1] .- by_eye

        a = a.^2

        a = a.^0.5

        optimum = argmin(a) + 1

        too_high_thresholds = hcat(repeat(["NA"],length(threshold_limit)),threshold_limit)

        entry_by_eye = hcat(by_eye,"by_eye")

        result = vcat(result, too_high_thresholds)

        result = vcat(result, entry_by_eye)

        writedlm(string(dir, "threshold_results.csv"),result,',',)

        println(" All tresholds evaluated.")

end


function particle_counter(
    file,
    image_format,
    intensity_filter::Float64,
    intensity_filter_mode,
    filter_color,
    vertical_steps,
    use_kernel,
    crop,
    x1,x2,y1,y2)

    if image_format == "JPG"

        results = image_load_and_filter_jpg(
            file,
            intensity_filter,
            intensity_filter_mode,
            filter_color,
            crop,
            x1,
            x2,
            y1,
            y2
        )
    end

    if image_format == "tif"

        results = image_load_and_filter_tif(
            file,
            intensity_filter,
            intensity_filter_mode,
            filter_color,
            use_kernel,
            crop,
            x1,
            x2,
            y1,
            y2
        )

    end



    img_bit = results[1]

    pic_height = results[2]

    pic_width = results[3]

    # subset only the white pixels, and store the addresses

    white_pixels = findall(img_bit .> 0)

    #  converts the image data, into a one dimensional address of all white pixels

    print("\n \n ...reducing adresses...")

    column_vector = address_dim_reduce(white_pixels, pic_height)

    # cuts columns into connected white_pixel stripes


    print("\n ...cutting into stripes...")


    stripes = cut_stripes(column_vector)

    print("\n ...connecting stripes to the right...")

    stripes =  v_step_connector(stripes,pic_height, vertical_steps)

    print("\n ...connecting stripes to the left...")

    particles = reverse_connecting_stripes(stripes, pic_height)


    return particles, pic_height, pic_width

end


function column_seperation(columns_wp, wp_convert, pic_height)

    column_vector = Vector{Vector{Int}}(undef, 0)

    # calculate one d address

    one_d = wp_convert[1, :] .+ (wp_convert[2, :] .- 1) .* pic_height

    counts = countmap(wp_convert[2, :])

    from = 1

    adder = 0

    for i in columns_wp::Array{Int64,1}

        # findall values that are in column i
        result = counts[i] + adder

        wp = one_d[from:result]

        wp = [wp]

        #print("\n wp ", typeof(wp))

        column_vector = vcat(column_vector, wp)

        from = result + 1

        adder = result

    end

    return column_vector

end

function address_dim_reduce(white_pixels, pic_height)


    # function that converts caretsianIndex into an array, than the column values are easily callable

    as_ints(a::AbstractArray{CartesianIndex{L}}) where {L} =
        reshape(reinterpret(Int, a), (L, size(a)...))

    wp_convert = as_ints(white_pixels)

    #print("\n format wp_convert", wp_convert[2, :])

    columns_wp = wp_convert[2, :]

    columns_wp = unique(columns_wp)

    #print("\n allocation vector", summary(column_vector))

    # here it checks for the white_pixels in which column they sit and converts them to the
    # 1d address of the position and devides them into sets of single column


    column_vector = column_seperation(columns_wp, wp_convert, pic_height)


    return column_vector

end

function cut_stripes(column_vector)

    stripes = Vector{Vector{Int}}(
        undef,
        length(column_vector) * maximum(length.(column_vector)),
    )

    # position in results vector

    position = 1

    for s = 1:length(column_vector)

        # number of white pixels in the called column

        max_row = length(column_vector[s])

        # loop for the column

        adder = 0

        for p = 1:length(column_vector[s])
            stripe = Vector{Int}(undef, length(column_vector[s]))

            p = p + adder

            if p > length(column_vector[s])

                break

            end

            stripe[1] = column_vector[s][p]

            n_p_adder = 1

            for n_p = (p+1):length(column_vector[s])

                if column_vector[s][p] + n_p_adder == column_vector[s][n_p]

                    stripe[1+n_p_adder] = column_vector[s][n_p]

                    adder = adder + 1

                    n_p_adder = n_p_adder + 1

                else

                    break

                end

            end

            # cut the vector to only def values

            stripe = stripe[1:n_p_adder]

            stripes[position] = stripe

            position = position + 1

        end

    end

    stripes = stripes[1:position-1]

    return stripes

end


# the following functions checks if stripes of white pixels are connected, if they are they get
# aggregated in one array

function connecting_stripes(stripes, pic_height)

    a = 1

    b = 0

    while a != b

        a = length(stripes)

        for i = 1:length(stripes)

            s = 1

            while i + s <= length(stripes)

                if stripes[i][end]::Int64 + pic_height::Int64 <
                           stripes[i+s][1]::Int64

                            break

                end

                y = intersect(stripes[i] .+ pic_height::Int64, stripes[i + s])

                if length(y) != 0

                    stripes[i] = vcat(stripes[i], stripes[i+s]::Array{Int64,1})

                    stripes = stripes[1:end .!= i+s, :]

                    s = s - 1

                end

                s = s + 1

            end

        end

        b = length(stripes)

    end

    return stripes

end

function connecting_stripes_slow(stripes, pic_height)

    # idea: maybe adjust this function to only check stripes for connections with length above X

    a = 1

    b = 0

    while a != b

        a = length(stripes)

        for i = 1:length(stripes)

            s = 1

            while i + s <= length(stripes)

                y = intersect(stripes[i] .+ pic_height::Int64, stripes[i + s])

                if length(y) != 0

                    stripes[i] = vcat(stripes[i], stripes[i+s]::Array{Int64,1})

                    stripes = stripes[1:end .!= i+s, :]

                    s = s - 1

                end

                s = s + 1

            end

        end

        b = length(stripes)

    end

    return stripes

end


function reverse_connecting_stripes(stripes, pic_height)

    a = 1

    b = 0

    while a != b

        a = length(stripes)

        for i = 1:length(stripes)

            s = 1

            while i + s <= length(stripes)


                if maximum(stripes[i])::Int64 + pic_height::Int64 <
                           stripes[i+s][1]::Int64

                            break

                end

                y = intersect(stripes[i] .- pic_height::Int64, stripes[i + s])

                if length(y) != 0

                    stripes[i] = vcat(stripes[i], stripes[i+s]::Array{Int64,1})

                    stripes = stripes[1:end .!= i+s, :]

                    s = s - 1

                end

                s = s + 1

            end

        end

        b = length(stripes)

    end

    return stripes

end

function v_step_connector(stripes, pic_height, steps)

    stepper = round(length(stripes)/steps, digits = 0)

    boundaries = convert.(Int,round.(collect(1:stepper:length(stripes)), digits = 0))

    step_test = Vector{Vector{Int}}(undef,0)

    for i in 1:length(boundaries)

        if i == length(boundaries)

        step_stripes = connecting_stripes(stripes[boundaries[i]+1:end], pic_height)

        step_test = vcat(step_test, step_stripes)

        break

        end

        if i == 1

        step_stripes = connecting_stripes(stripes[boundaries[1]:boundaries[i+1]], pic_height)

        step_test = vcat(step_test, step_stripes)

        end

        if i > 1

        step_stripes = connecting_stripes(stripes[boundaries[i]+1:boundaries[i+1]], pic_height)

        step_test = vcat(step_test, step_stripes)

        end

    end

    particles = connecting_stripes(step_test, pic_height)

    return particles

end

function h_step_connector(stripes, pic_height,pic_width, steps)

    positions = particle_positions(pic_height, pic_width, stripes)

    stepper = pic_height/steps

    boundaries =
        convert.(Int, round.(collect(1:stepper:pic_height), digits = 0))


    step_test = Vector{Vector{Int}}(undef, 0)

    # this for loop uses the boundaries to call connecting stripes on parts of the whole package of stripes

    for i = 1:length(boundaries)


        if i == length(boundaries)

            h_step = findall(x -> x >= boundaries[end]+1 && x <= pic_height, positions[2])

            #println(" from ", boundaries[end]+1, " to ", pic_height)

            step_stripes = connecting_stripes(stripes[h_step], pic_height)

            step_test = vcat(step_test, step_stripes)

            break

        end

        if i == 1

            h_step = findall(x -> x >= boundaries[1] && x <= boundaries[i + 1], positions[2])

            #println(" from ", boundaries[1], " to ", boundaries[i + 1])

            step_stripes = connecting_stripes(
                stripes[h_step],
                pic_height,
            )

            step_test = vcat(step_test, step_stripes)

        end

        if i > 1

            h_step = findall(x -> x >= boundaries[i]+1 && x <= boundaries[i + 1], positions[2])

            #println(" from ", boundaries[i] + 1, " to ", boundaries[i + 1])

            step_stripes = connecting_stripes(
                stripes[h_step],
                pic_height,
            )

            step_test = vcat(step_test, step_stripes)

        end

    end

    return step_test

end


function image_load_and_filter_jpg(
    file,
    intensity_filter,
    intensity_filter_mode,
    filter_color,
    crop,
    x1,
    x2,
    y1,
    y2)

    img = load(file)

    # convert the jpg into channels/arrays

    img_array = channelview(img)


    if filter_color == "green"

        # subset only the dimension containing the green values

        img_monochrome = img_array[2, :, :]

        # convert the values from UInt8 to Float64 or something you can deal with

        img_monochrome = convert(Array{Float64}, img_monochrome)

        if crop == true

            img_monochrome = img_monochrome[y1:y2,x1:x2]

            println("\n Only analysed cropped image!")

        end

        a = minimum(img_monochrome)

        b = round(maximum(img_monochrome), digits = 2)

        print("\n // Max green value:", b, "//")

    end

    if filter_color == "red"

        # subset only the dimension containing the red values

        img_monochrome = img_array[1, :, :]

        img_monochrome = convert(Array{Float64}, img_monochrome)

        if crop == true

            img_monochrome = img_monochrome[y1:y2,x1:x2]

        end

        a = minimum(img_monochrome)

        b = round(maximum(img_monochrome), digits = 2)

        print("\n // Max red value:", b, "//")

    end

    if filter_color == "blue"

        # subset only the dimension containing the red values

        img_monochrome = img_array[3, :, :]

        img_monochrome = convert(Array{Float64}, img_monochrome)

        if crop == true

            img_monochrome = img_monochrome[y1:y2,x1:x2]

        end

        a = minimum(img_monochrome)

        b = round(maximum(img_monochrome), digits = 2)

        print("\n // Max blue value:", b, "//")

    end



    # filter the image noise by a treshold; the array is converted into a bit array, that's perfect:)!

    if intensity_filter_mode == "relative"

        img_bit = img_monochrome .> b * intensity_filter

        print("\n // Threshold for this image:", b * intensity_filter, "//")

    end

    if intensity_filter_mode == "absolute"

        img_bit = img_monochrome .> intensity_filter

        print("\n // Threshold for this image:", intensity_filter, "//")

    end

    pic_height = length(view(img_bit, :, 1))

    pic_width = length(view(img_bit, 1, :))

    return img_bit, pic_height, pic_width

end

function image_load_and_filter_tif(
    file,
    intensity_filter,
    intensity_filter_mode,
    filter_color,
    use_kernel,
    crop,
    x1,
    x2,
    y1,
    y2)

    img = load(file)

    # convert the jpg into channels/arrays

    img_array = channelview(img)


    if filter_color == "green"

        # subset only the dimension containing the green values

        img_monochrome = img_array[2, :, :,1]

        # convert the values from UInt8 to Float64 or something you can deal with

        img_monochrome = convert(Array{Float64}, img_monochrome)

        if crop == true

            img_monochrome = img_monochrome[y1:y2,x1:x2]

            println("\n Only analysed cropped image!")

        end

        if use_kernel == true

            a = repeat([1/9], 9)

            a = reshape(a,3,3)

            kernel = OffsetArray(a,-1:1,-1:1)

            img_monochrome = imfilter(img_monochrome, kernel)

            println("\n Used an 3x3 averaging kernel!")

        end

        a = minimum(img_monochrome)

        b = round(maximum(img_monochrome), digits = 2)

        print("\n // Max green value:", b, "//")

    end

    if filter_color == "red"

        # subset only the dimension containing the red values

        img_monochrome = img_array[1, :, :,1]

        img_monochrome = convert(Array{Float64}, img_monochrome)

        if crop == true

            img_monochrome = img_monochrome[y1:y2,x1:x2]

        end

        if use_kernel == true

            a = repeat([1/9], 9)

            a = reshape(a,3,3)

            kernel = OffsetArray(a,-1:1,-1:1)

            img_monochrome = imfilter(img_monochrome, kernel)

            println("\n Used an averaging kernel!")

        end

        a = minimum(img_monochrome)

        b = round(maximum(img_monochrome), digits = 2)

        print("\n // Max red value:", b, "//")

    end

    if filter_color == "blue"

        # subset only the dimension containing the red values

        img_monochrome = img_array[3, :, :,1]

        img_monochrome = convert(Array{Float64}, img_monochrome)

        if crop == true

            img_monochrome = img_monochrome[y1:y2,x1:x2]

        end

        if use_kernel == true

            a = repeat([1/9], 9)

            a = reshape(a,3,3)

            kernel = OffsetArray(a,-1:1,-1:1)

            img_monochrome = imfilter(img_monochrome, kernel)

            println("\n Used an averaging kernel!")

        end

        a = minimum(img_monochrome)

        b = round(maximum(img_monochrome), digits = 2)

        print("\n // Max blue value:", b, "//")

    end



    # filter the image noise by a treshold; the array is converted into a bit array, that's perfect:)!

    if intensity_filter_mode == "relative"

        img_bit = img_monochrome .> b * intensity_filter

        print("\n // Threshold for this image:", b * intensity_filter, "//")

    end

    if intensity_filter_mode == "absolute"

        img_bit = img_monochrome .> intensity_filter

        print("\n // Threshold for this image:", intensity_filter, "//")

    end

    pic_height = length(view(img_bit, :, 1))

    pic_width = length(view(img_bit, 1, :))

    return img_bit, pic_height, pic_width

end

function original_crop_png(file,x1, x2, y1, y2)

    img = load(file)

    img_crop = img[y1:y2,x1:x2,1]

    save(string(file[begin:end-4],"_img_crop.png"), img_crop)

end


function re_image_color(width, height, particles)

    img = zeros(width * height) # creates an "empty" bit-array

    for s in 1:length(particles)

        color = rand(0.3:0.01:1.0,1)

        for i in particles[s]

            img[i] = color[1]

        end

    end

    img_post = reshape(img, height, :)

    return colorview(Gray, img_post)

end

function re_image(width, height, particles)

    b = vcat(particles...) #  a vector with every white-pixel in a 1d-address format

    img = falses(width * height) # creates an "empty" bit-array

    for i in b
        img[i] = 1

    end

    img_post = reshape(img, height, :)

    return colorview(Gray, img_post)

end


function particle_positions(pic_height, pic_width, particles)

    points = minimum.(particles)

    x_positions = cld.(points, pic_height)

    y_positions = mod1.(points, pic_height)

    return x_positions, y_positions

end

# extracts the metadata of an image via shellout using exiftool

function get_metadata_timestamp(img, image_format)

    if image_format == "JPG"

        mdata = read(`exiftool -SonyDateTime $img`, String)

    end

    if image_format == "tif"

        mdata = read(`exiftool -createdate $img`, String)

    end

    timestamp = mdata[end-19:end-1]

    return timestamp

end

function get_metadata_shutter_speed(img, image_format)

    if image_format == "JPG"

        mdata = read(`exiftool -shutterspeed $img`, String)

    end

    if image_format == "tif"

        mdata = read(`exiftool -exposuretime $img`, String)

    end

    shutter_speed = mdata[35:end-1]

    return shutter_speed

end

function get_metadata_aperture(img)

    mdata = read(`exiftool -aperture $img`, String)

    aperture = mdata[35:end-1]

    return aperture

end

function get_metadata_ISO(img, image_format)

    if image_format == "JPG"

        mdata = read(`exiftool -ISO $img`, String)

        ISO = mdata[35:end-1]

    end

    if image_format == "tif"

        mdata = read(`exiftool -recommendedexposureindex $img`, String)

        ISO = mdata[35:end-1]

    end

    return ISO

end

function get_metadata_resolution(img, image_format)

    if image_format == "JPG"

        mdata = read(`exiftool -ImageSize $img`, String)

        width = mdata[end-9:end-6]

        height = mdata[end-4:end-1]

    end

    if image_format == "tif"

        mdata = read(`exiftool -imageheight $img`, String)

        height = mdata[35:end-1]

        mdata = read(`exiftool -imagewidth $img`, String)

        width = mdata[35:end-1]

    end

    return width, height

end

function show_evaluated_particles(
    file,
    pic_width,
    pic_height,
    filtered_particles,
    mark_size,
    intensity_filter)

    img = load(file)

    # convert the jpg into channels/arrays

    img_array = channelview(img)

    # subset only the dimension containing the green values

    img_monochrome = img_array[2, :, :,1]

    # convert the values from UInt8 to Float64 or something you can deal with

    img_monochrome = convert(Array{Float64}, img_monochrome)

    marking_pixels =
        calculate_marker(filtered_particles, pic_width, pic_height, mark_size)

    print("\n // start marking //")

    marked_img = mark_particles(img_monochrome, marking_pixels, pic_height)

    print("\n // built marked image //")

    save(
        string(file[1:end-4], "_", intensity_filter, "_marked.png"),
        marked_img,
    )

end

function mark_particles(img, marker, pic_height)

    b = vcat(img...) #  a vector with every white-pixel in a 1d-address format


    for i in marker
        b[i] = 1

    end

    print("\n reshape ")

    img_post = reshape(b, pic_height, :)

    print("\n restructure ")

    marked_img = colorview(Gray, img_post)

    return marked_img

end

function calculate_marker( filtered_particles,
    pic_width,
    pic_height,
    mark_size,)

    #center_points = round.(median.(filtered_particles), digits = 0)

    center_points = round.(minimum.(filtered_particles), digits = 0)

    mark_pixels = Int64[]

    for c in center_points
        start_point = c - mark_size / 2 - mark_size / 2 * pic_height

        mark_pixels = vcat(mark_pixels, start_point)

        for i = 1:mark_size-1
            start = start_point

            start = start + i

            mark_pixels = vcat(mark_pixels, start)

        end

        for i = 1:mark_size-1
            start = start_point

            start = start + i + pic_height * mark_size

            mark_pixels = vcat(mark_pixels, start)

        end

        for i = 1:mark_size-1
            start = start_point

            start = start + pic_height * i

            mark_pixels = vcat(mark_pixels, start)

        end

        for i = 1:mark_size-1
            start = start_point

            start = start + pic_height * i + mark_size

            mark_pixels = vcat(mark_pixels, start)

        end

    end

    mark_pixels = floor.(Int, mark_pixels)

    mark_pixels = filter((x) -> x > 0, mark_pixels)

    mark_pixels = filter((x) -> x < pic_width * pic_height, mark_pixels)

    return mark_pixels

end


function perimeter_connection(particles, pic_height)

    pixel_perimeters = Vector{Vector{Int}}(undef, 0)

    a = Int64[]

    b = Int64[]

    c = Int64[]

    d = Int64[]

    for particle in particles
        pixel_perimeter = Int64[]

        for pixel in particle
            a = findall(x -> x == pixel + 1, particle)

            b = findall(x -> x == pixel - 1, particle)

            c = findall(x -> x == pixel + pic_height, particle)

            d = findall(x -> x == pixel - pic_height, particle)

            result = vcat(length(a), length(b), length(c), length(d))

            pixel_perimeter = vcat(pixel_perimeter, result)

        end

        pixel_perimeters = vcat(pixel_perimeters, [pixel_perimeter])

    end

    return pixel_perimeters

end

function aspect_ratio(particles, pic_height)

    # pre allocating two_D possible

    aspect_ratio_results_vh = Vector{Float64}(undef, length(particles))

    aspect_ratio_results_diagonal = Vector{Float64}(undef, length(particles))

    extremes = Vector{Float64}(undef, 4)

    counter = 1

    for particle in particles
        distance_BR = Vector{Float64}(undef, length(particle))

        distance_TR = Vector{Float64}(undef, length(particle))

        distance_BL = Vector{Float64}(undef, length(particle))

        distance_TL = Vector{Float64}(undef, length(particle))

        extremes[1] = minimum(cld.(particle, pic_height)) # left

        extremes[2] = minimum(mod1.(particle, pic_height)) # top

        extremes[3] = maximum(cld.(particle, pic_height)) # right

        extremes[4] = maximum(mod1.(particle, pic_height)) # bottom

        x = cld.(particle, pic_height)

        y = mod1.(particle, pic_height)

        two_D = hcat.(x, y)

        # calculate distance to corners

        for i = 1:length(two_D)
            distance_BR[i] = sqrt(
                (two_D[i][1] - extremes[3])^2 + (two_D[i][2] - extremes[4])^2,
            )

        end

        for i = 1:length(two_D)
            distance_TR[i] = sqrt(
                (two_D[i][1] - extremes[3])^2 + (two_D[i][2] - extremes[2])^2,
            )

        end

        for i = 1:length(two_D)
            distance_BL[i] = sqrt(
                (two_D[i][1] - extremes[1])^2 + (two_D[i][2] - extremes[4])^2,
            )

        end

        for i = 1:length(two_D)
            distance_TL[i] = sqrt(
                (two_D[i][1] - extremes[1])^2 + (two_D[i][2] - extremes[2])^2,
            )

        end

        vertical_distance = extremes[4] - extremes[2]

        horizontal_distance = extremes[3] - extremes[1]

        TL = two_D[findall(x -> x == minimum(distance_TL), distance_TL)][1]

        TR = two_D[findall(x -> x == minimum(distance_TR), distance_TR)][1]

        BL = two_D[findall(x -> x == minimum(distance_BL), distance_BL)][1]

        BR = two_D[findall(x -> x == minimum(distance_BR), distance_BR)][1]

        diagonal_distance_TL_BR =
            round(sqrt((TL[1] - BR[1])^2 + (TL[2] - BR[2])^2), digits = 2)

        diagonal_distance_BL_TR =
            round(sqrt((BL[1] - TR[1])^2 + (BL[2] - TR[2])^2), digits = 2)

        aspect_ratio_vh =
            round(vertical_distance / horizontal_distance, digits = 1)

        aspect_ratio_diagonal =
            round(diagonal_distance_BL_TR / diagonal_distance_TL_BR, digits = 1)

        aspect_ratio_results_vh[counter] = aspect_ratio_vh

        aspect_ratio_results_diagonal[counter] = aspect_ratio_diagonal

        counter = counter + 1

    end

    results = hcat(aspect_ratio_results_vh, aspect_ratio_results_diagonal)

    return results

end

function write_hist(table_hist, particles, postfix, timestamp)

    sizes = length.(particles)

    if length(sizes) == 0

        hist_data = hcat(NA, timestamp, postfix)

    else

        hist_timestamp = repeat([timestamp], length(sizes))

        hist_type = repeat([postfix], length(sizes))

        table_entry = hcat(sizes, hist_timestamp, hist_type)

    end

    # returns just the basic hist entry

    return table_entry

end

function particle_data(
    dir,
    intensity_filter;
    filter_color::String = "green",
    intensity_filter_mode::String = "absolute",
    vertical_steps::Int64 = 10,
    filter_image::Bool = false,
    timestamp_mode::String = "metadata",
    mark_mode::Bool = false,
    mark_size::Int64 = 40,
    from::Int64 = 1,
    to::Int64 = 0,
    circularity_analysis::Bool = false,
    image_format::String = "tif",
    aspect_ratio_mode::Bool = false,
    position_calculation::Bool = false,
    dont_write_results::Bool = false,
    use_kernel::Bool = false,
    crop::Bool = false,
    x1::Int64 = 1,
    x2::Int64 = 100,
    y1::Int64 = 1,
    y2::Int64 = 100,
)

    files = readdir(dir, join = true)

    if image_format == "JPG"

        filter!(s -> occursin(r".JPG", s), files)

    end

    if image_format == "tif"

        filter!(s -> occursin(r".tif", s), files)

    end

    print(
        "\n // I found ",
        length(files),
        " files",
        "\n first file: \n ",
        files[1],
        " last file: \n",
        files[end],
        " //",
    )

    if timestamp_mode == "metadata"

        print(
            "\n // Time information will be extracted from image metadata  // ",
        )

        table_hist = hcat("size", "time", "type")

        table_positions = hcat("X_positions", "y_positions", "time")

        if circularity_analysis == true

            table_hist = hcat(table_hist, "circularity")

        end

        if aspect_ratio_mode == true

            table_hist = hcat(table_hist, "ar_vh", "ar_dia")

        end

    end

    ###############################
    # batch process starts here #
    ###############################


    if to > 0

        to = abs(to - length(files))

    end

    progress = 1

    for file in files[from:end-to]
        print("\n \n start processing: ", file)

        # extract the timestamp from metadata

        if timestamp_mode == "metadata"

            timestamp = get_metadata_timestamp(file, image_format)

        end

        results = particle_counter(
            file,
            image_format,
            intensity_filter,
            intensity_filter_mode,
            filter_color,
            vertical_steps,
            use_kernel,
            crop,
            x1,
            x2,
            y1,
            y2,
        )

        # results are saved from the particle_counter function for use in other functions

        particles = results[1]

        # height of the pic

        height = results[2]

        # width of the pic

        width = results[3]

        return particles

    end



end
