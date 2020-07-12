configfile: "config_germany.yaml"

subworkflow pypsade:
    workdir: "pypsa-eur"
    snakefile: "pypsa-eur/Snakefile"
    configfile: "config_germany.yaml"


rule sanitize_network:
    input:
        network = pypsade('results/networks/elec_s_{clusters}_ec_lvopt_Ep.nc')
    output: 
        network = 'resources/german_network_{clusters}.nc'

rule plot_german_network:
    input:
        network = 'resources/german_network_{clusters}.nc',
        regions = pypsade('resources/regions_onshore_elec_s_{clusters}.geojson')
    output: 'figures/german_network_{clusters}.png'
    script: 'code/plot_network.py'


rule allocate_network:
    input:
        network = 'resources/german_network_{clusters}.nc'
    output:
        costs = 'resources/allocated_costs_{clusters}_{method}_{power}.nc',
        payments = 'resources/allocated_payments_{clusters}_{method}_{power}.nc'
    script: 'code/allocate_network.py'


rule plot_total_costs:
    input:
        costs = 'resources/allocated_costs_{clusters}_ptpf_net.nc'
    output: 'figures/total_costs_{clusters}.png'
    script: 'code/plot_total_costs.py'

rule plot_price_maps:
    input:
        network = 'resources/german_network_{clusters}.nc',
        regions = pypsade('resources/regions_onshore_elec_s_{clusters}.geojson'),
        payments = 'resources/allocated_payments_{clusters}_{method}_{power}.nc'
    output:
        folder = directory('figures/price_maps_{clusters}_{method}_{power}')
    script: 'code/plot_price_map.py'



# Local Variables:
# mode: python
# End:
