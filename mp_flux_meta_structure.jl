function mp_flux(
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
    write_results::Bool = true,
    use_kernel::Bool = false,
    crop::Bool = false,
    x1::Int64 = 1,
    x2::Int64 = 100,
    y1::Int64 = 1,
    y2::Int64 = 100,
    area_density::Bool = false,
    N_f::Int64 = 1,
    kernel_length::Int64 = 285,
    center_of_mass_image::Bool = false
)

    print("\n Welcome to mp_flux. \n \n")

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

        print("\n // Time information will be extracted from image metadata  // ")

        table_hist = hcat("size", "time", "type")

        table_positions = hcat("X_positions", "y_positions", "time")

        if circularity_analysis == true

            table_hist = hcat(table_hist, "circularity")

        end

        if aspect_ratio_mode == true

            table_hist = hcat(table_hist, "ar_vh", "ar_dia")

        end

        if area_density == true

            area_density_table = hcat("time","mean","max","N")

        end

    end

    if timestamp_mode == "postfix"

        table_hist = hcat("size", "postfix", "type")

        table_positions = hcat("X_positions", "y_positions", "postfix")

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

        # extract the timestamp from filename

        if timestamp_mode == "postfix"

            timestamp = split(file, "-")[end]

            timestamp = split(timestamp, ".")[1]

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
            x1,x2,y1,y2
        )

        # results are saved from the particle_counter function for use in other functions

        particles = results[1]

        # height of the pic

        pic_height = results[2]

        # width of the pic

        pic_width = results[3]

        basic_hist_entry =
            write_hist(table_hist, particles, "pre", timestamp)

        if circularity_analysis == true

            print("\n CIRCULARITY analysis...")

            a = perimeter_connection(particles, height)

            perimeters = length.(a) - sum.(a)

            sizes = length.(particles)

            circularity =
                round.(sqrt.((4 * pi .* sizes) ./ perimeters .^ 2), digits = 2)

            # adds the information to the entry

            basic_hist_entry = hcat(basic_hist_entry, circularity)

        end

        if aspect_ratio_mode == true

            print("\n start calculating aspect ratios: ")

            aspect_ratios = aspect_ratio(particles, height)

            basic_hist_entry = hcat(basic_hist_entry, aspect_ratios)

        end

        # adding the gathered information to the table

        table_hist = vcat(table_hist, basic_hist_entry)

        if position_calculation == true

            positions = particle_positions(height, pic_width, particles)

            pos_timestamp = repeat([timestamp], length(positions[1]))

            positions_entry = hcat(positions[1], positions[2], pos_timestamp)

            table_positions = vcat(table_positions, positions_entry)

        end


        # write hist data into table after the filter is applied

        print("\n Processed: ", file)

        # re-image settings

        if filter_image == true

            img = re_image_color(pic_width, pic_height, particles)

            save(string(file[1:end-4], "_", intensity_filter, ".png"), img)

            print("\n // image of threshold surviving pixels saved //")

        end

        if mark_mode == true

            print("\n ... building a marked version of the original ...")

            show_evaluated_particles(
                file,
                pic_width,
                pic_height,
                particles,
                mark_size,
                intensity_filter,
            )

        end

        if area_density == true

            println("... calculating center of mass ...")

        center_results =  center_of_mass_address(particles, N_f, pic_height)

        writedlm(string(file[1:end-4], "_", "area_density", ".csv"), center_results,',')

        end

        if center_of_mass_image == true

            re_image_center_of_mass(pic_width,pic_height,particles,ad_results[3],dir, file)

            println("... saved a center of mass image ...")

        end

        print("\n (", progress, "/", length(files[from:end-to]), ")")

        progress = progress + 1

    end

    # create an info table

    header = hcat("type", "value")

    row_1 = hcat("intenstiy_filter", intensity_filter)

    res = get_metadata_resolution(files[end], image_format)

    row_2 = hcat("width", res[1])

    row_3 = hcat("height", res[2])

    row_4 = hcat("ISO", get_metadata_ISO(files[end], image_format))

    row_5 = hcat("aperture", get_metadata_aperture(files[end]))

    row_6 = hcat("shutter_speed", get_metadata_shutter_speed(files[end], image_format))

    info_table = vcat(header, row_1,row_2, row_3, row_4, row_5, row_6)

    ##############################################
    # save tables at the directory of the images #
    ##############################################

    if write_results == true

        if crop == true

            writedlm(
                string(dir, "/","cropped_", "particle_hist_data", "_", intensity_filter, ".csv"),
                table_hist,
                ',',
            )

            writedlm(
                string(dir, "/","cropped_", "setup_data_info", "_", intensity_filter, ".csv"),
                info_table,
                ',',
            )

            if position_calculation == true

                writedlm(
                    string(dir, "/","cropped_", "positions", "_", intensity_filter, ".csv"),
                    table_positions,
                    ',',
                )

            end

        else

            writedlm(
                string(dir, "/", "particle_hist_data", "_", intensity_filter, ".csv"),
                table_hist,
                ',',
            )

            writedlm(
                string(dir, "/", "setup_data_info", "_", intensity_filter, ".csv"),
                info_table,
                ',',
            )

            if position_calculation == true

                writedlm(
                    string(dir, "/", "positions", "_", intensity_filter, ".csv"),
                    table_positions,
                    ',',
                )

            end

            if area_density == true

                writedlm(
                    string(dir, "/", "area_density.csv"),
                    area_density_table,
                    ',',
                )

            end

        end

        print("\n ##### FINISHED ALL FILES #####")

        print("\n Results are saved in ", dir)

    else

        println("\n DID NOT SAVE RESULTS, DUE TO OPTION write_results")

    end

end
