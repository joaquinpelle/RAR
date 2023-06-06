using Skylight
using CairoMakie
using Printf
using DataInterpolations

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

inner_radii = Dict("r0"=>1.654835e+06, "rsat"=>1.980475e+12)

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

        image_plane = ImagePlane(distance = 1e-5*r_sup,
                                observer_inclination_in_degrees = ξ,
                                horizontal_side = srsat*rsat,
                                vertical_side = srsat*rsat,
                                horizontal_number_of_pixels = Nres,
                                vertical_number_of_pixels = Nres)

        model = Skylight.ShakuraSunyaevDisk(inner_radius=rd_in, outer_radius=rd_out, alpha=0.5, M1=M1)
                
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                                camera = image_plane,
                                                observed_times = [0.0],
                                                radiative_model = model,
                                                unit_mass_in_solar_masses=model.M1)
        initial_data = load_initial_data_from_hdf5("io/schw/$(filename).h5")
        output_data = load_output_data_from_hdf5("io/schw/$(filename).h5", 1)
        # configurations = load_configurations_from_hdf5("io/schw/$(filename).h5")            
        
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
        save("plots/schw/line_broadening/$(filename)_bins$num_bins.png", fig; dpi=300)
    end
end