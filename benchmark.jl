using Skylight
using BenchmarkTools

M1 = 1.02e7 #Unit mass in solar masses
energies = [200] #energies in keV of the fermion
inclinations = [42] #inclinations in degrees of the observer
Nres = 100 #number of pixels per side of the image
srsat = 75 #side of the image plane in units of Rsat

rsat = 1.980475e+12
rsat = CGS_to_geometrized(rsat, Dimensions.length, M1=M1) 
rd_out = 1e3*rsat

for E in energies

    spacetime = RARSpacetime(data_dir = "io/$(E)keV")
    r_inf, r_sup = radial_bounds(spacetime)

    rd_in = r_inf
    for ξ in inclinations
        
        camera = ImagePlane(distance = 1e-5*r_sup,
                                observer_inclination_in_degrees = ξ,
                                observation_times = [0.0],
                                horizontal_side = srsat*rsat,
                                vertical_side = srsat*rsat,
                                horizontal_number_of_pixels = Nres,
                                vertical_number_of_pixels = Nres)

        model = RARDisk(inner_radius=rd_in, outer_radius=rd_out, alpha=0.5, M1=M1)
                
        configurations = VacuumOTEConfigurations(spacetime=spacetime,
                                                camera = camera,
                                                radiative_model = model,
                                                unit_mass_in_solar_masses=model.M1)

        initial_data = initialize(configurations)

        cb, cb_params = callback_setup(configurations) #... or, define your own cb and cb_params

        integrate(initial_data, configurations, cb, cb_params; method=VCABM(), reltol=1e-13, abstol=1e-21)

        
    end
end
