using Skylight
using CairoMakie
using Printf

M1 = 1.02e7 #Unit mass in solar masses
energies = [200] #energies in keV of the fermion
inclinations = [5, 25, 45, 65, 85] #inclinations in degrees of the observer
Nres = 1000 #number of pixels per side of the image
srsat = 75 #side of the image plane in units of Rsat

inner_radii = Dict("r0"=>1.654835e+06, "rsat"=>1.980475e+12)

function myprofile(position, spacetime, model)
    r = radius(position, spacetime)
    return 1/r^2
end

num_bins = 40

for (rname, rd_in) in inner_radii
    inner_radii[rname] = CGS_to_geometrized(rd_in, Dimensions.length, M1=M1)
end
    
rsat = inner_radii["rsat"]
rd_out = 1e3*rsat

spacetime = SchwarzschildSpacetimeSphericalCoordinates(M=1.0)
rd_in = isco_radius(spacetime)

for E in energies

    rar_spacetime = RARSpacetime(data_dir = "io/$(E)keV")
    r_inf, r_sup = radial_bounds(rar_spacetime)

    for ξ in inclinations
            
        ξstr = string(@sprintf("%02d", ξ)) 
        Nstr = string(@sprintf("%04d", Nres))
        filename = "$(E)keV_ideg$(ξstr)_s$(srsat)rsat_N$(Nstr)"

        camera = ImagePlane(distance = 1e-5*r_sup,
                                observer_inclination_in_degrees = ξ,
                                observation_times = [0.0],
                                horizontal_side = srsat*rsat,
                                vertical_side = srsat*rsat,
                                horizontal_number_of_pixels = Nres,
                                vertical_number_of_pixels = Nres)

        model = Skylight.ShakuraSunyaevDisk(inner_radius=rd_in, outer_radius=rd_out, alpha=0.5, M1=M1)
                
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                                camera = camera,
                                                radiative_model = model,
                                                unit_mass_in_solar_masses=model.M1)

        initial_data = initialize(configurations)

        cb, cb_params = callback_setup(configurations, rhorizon_bound=2e-1) #... or, define your own cb and cb_params

        run = integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)
        output_data = run.output_data

	    save_to_hdf5("io/schw/$(filename).h5", configurations, initial_data, run)        
        
        #Bolometric intensity image
        Iobs, q = observed_bolometric_intensities(initial_data, output_data, configurations)

        xs,ys = axes_ranges(configurations.camera)

        zs = grid_view(Iobs, configurations)

        fig = Figure(font = "Times New Roman")
        ax = Axis(fig[1,1], xlabel=L"\alpha/(GM/c^2)", ylabel=L"\beta/(GM/c^2)", ylabelsize = 26, xlabelsize = 26) 
        hmap = heatmap!(xs, ys, zs/maximum(zs); colormap=:gist_heat, interpolate=true)
        Colorbar(fig[:, end+1], hmap, label=L"I \text{(arbitrary)}", labelsize=26, width = 15, ticksize = 18, tickalign = 1)
        colsize!(fig.layout, 1, Aspect(1, 1.0))
        colgap!(fig.layout, 7)
        CairoMakie.save("plots/schw/bolometric/$(filename).png", fig)
        
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
        CairoMakie.save("plots/schw/specific/$(filename).png", fig)

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
        CairoMakie.save("plots/schw/line_broadening/$(filename)_bins$num_bins.png", fig; dpi=300)
    end
end
