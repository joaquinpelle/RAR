using Skylight
using CairoMakie
using Printf

M1 = 1.02e7 #Unit mass in solar masses
energies = [200] #energies in keV of the fermion
inclinations = [5, 25, 45, 65, 85] #inclinations in degrees of the observer
Nres = 1000 #number of pixels per side of the image
srsat = 75 #side of the image plane in units of Rsat

function myprofile(position, spacetime, model)
    r = radius(position, spacetime)
    return 1/r^2
end

num_bins = 40

rsat = 1.980475e+12
rsat = CGS_to_geometrized(rsat, Dimensions.length, M1=M1) 
rd_out = 1e3*rsat

for E in energies

    spacetime = RARSpacetime(data_dir = "io/$(E)keV")
    r_inf, r_sup = radial_bounds(spacetime)
    inner_radii = Dict("r0"=>r_inf, "rsat"=>rsat)

    for ξ in inclinations
        for (rname, rd_in) in inner_radii
            
            ξstr = string(@sprintf("%02d", ξ)) 
            Nstr = string(@sprintf("%04d", Nres))
            filename = "$(E)keV_ideg$(ξstr)_$(rname)_s$(srsat)rsat_N$(Nstr)"

            image_plane = ImagePlane(distance = 1e-5*r_sup,
                                    observer_inclination_in_degrees = ξ,
                                    horizontal_side = srsat*rsat,
                                    vertical_side = srsat*rsat,
                                    horizontal_number_of_pixels = Nres,
                                    vertical_number_of_pixels = Nres)

            model = RARDisk(inner_radius=rd_in, outer_radius=rd_out, alpha=0.5, M1=M1)
                    
            configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                                    image_plane = image_plane,
                                                    observed_times = [0.0],
                                                    radiative_model = model,
                                                    unit_mass_in_solar_masses=model.M1)

            initial_data = initialize(configurations)

            cb, cb_params = callback_setup(configurations) #... or, define your own cb and cb_params

            run = integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)

            save_to_hdf5("io/RAR/$(filename).h5", configurations, initial_data, run; mode="w")        
            
            output_data = run.output_data
            
            #Bolometric intensity image
            Iobs, q = observed_bolometric_intensities(initial_data, output_data, configurations)

            xs,ys = axes_ranges(camera)

            zs = grid_view(Iobs, configurations)

            fig = Figure(font = "Times New Roman")
            ax = Axis(fig[1,1], xlabel=L"\alpha/(GM/c^2)", ylabel=L"\beta/(GM/c^2)", ylabelsize = 26, xlabelsize = 26) 
            hmap = heatmap!(xs, ys, zs/maximum(zs); colormap=:gist_heat, interpolate=true)
            Colorbar(fig[:, end+1], hmap, label=L"I \text{(arbitrary)}", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
            colsize!(fig.layout, 1, Aspect(1, 1.0))
            colgap!(fig.layout, 7)
            CairoMakie.save("plots/RAR/bolometric/$(filename).png", fig)
            
            #Specific intensity image
            λ_EHT_Apr17 = 0.13
            ε = PhysicalConstants.h*PhysicalConstants.c/λ_EHT_Apr17

            Iobs, q = observed_specific_intensities(initial_data, output_data, configurations, [ε])

            zs = grid_view(Iobs, configurations; energy_index=1)

            fig = Figure(font = "Times New Roman")
            ax = Axis(fig[1,1], xlabel=L"\alpha/(GM/c^2)", ylabel=L"\beta/(GM/c^2)", ylabelsize = 26, xlabelsize = 26) 
            hmap = heatmap!(xs, ys, zs/maximum(zs); colormap=:gist_heat, interpolate=true)
            Colorbar(fig[:, end+1], hmap, label=L"I \text{(arbitrary)}", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
            colsize!(fig.layout, 1, Aspect(1, 1.0))
            colgap!(fig.layout, 7)
            CairoMakie.save("plots/RAR/specific/$(filename).png", fig)
            
            #Line emission spectrum

            binned_fluxes, bins = line_emission_spectrum(initial_data, output_data, configurations; emission_profile = myprofile, num_bins = num_bins)
            set_theme!(; fonts = (; regular = "Times New Roman"))

            # We calculate midpoints of x to use as x coordinates for y
            max_flux = maximum(binned_fluxes)
            bins_midpoints = 0.5*(bins[1:end-1] + bins[2:end])
            fig = Figure(resolution = (600, 400))
            ax = Axis(fig[1, 1], xlabel = L"E/E_0", ylabel = "Flux (arbitrary)", title = "Relativistic line broadening", titlefont=:regular)
            skl = lines!(ax, bins_midpoints, binned_fluxes/max_flux, linewidth = 3, color = :black)
            
            ax.titlesize = 22
            ax.xlabelsize = 22
            ax.ylabelsize = 22
            ax.xticklabelsize = 15
            ax.yticklabelsize = 15
            
            # Save the figure
            CairoMakie.save("plots/RAR/line_broadening/$(filename)_bins$num_bins.png", fig; dpi=300)

        end
    end
end
