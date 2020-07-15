configfile: "config_germany.yaml"


wildcard_constraints:
    clusters="[0-9]+m?|all",
    field="(bf|gf)"


subworkflow pypsade:
    workdir: "pypsa-eur"
    snakefile: "pypsa-eur/Snakefile"
    configfile: "config_germany.yaml"


rule solve_and_sanitize_german_network:
    input:
        network = pypsade('networks/elec_s_{clusters}_ec_lvopt_Ep.nc'),
        regions = pypsade('resources/regions_onshore_elec_s_{clusters}.geojson')
    output: 
        network = 'resources/de{clusters}{field}.nc',
        regions = 'resources/de{clusters}{field}_regions.geojson'
    script: 'code/solve_and_sanitize_network.py'


rule plot_network:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson'
    output: 'figures/{nname}.png'
    script: 'code/plot_network.py'


rule allocate_network:
    input:
        network = 'resources/{nname}.nc'
    output:
        costs = 'resources/costs_{nname}_{method}_{power}.nc',
        payments = 'resources/payments_{nname}_{method}_{power}.nc'
    script: 'code/allocate_network.py'


rule plot_total_costs:
    input:
        costs = 'resources/costs_{nname}_ptpf_net.nc'
    output: 'figures/total_costs_{nname}.png'
    script: 'code/plot_total_costs.py'

rule plot_price_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        payments = 'resources/payments_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/price_maps_{nname}_{method}_{power}')
    script: 'code/plot_price_map.py'

rule plot_cost_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        payments = 'resources/payments_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/cost_maps_{nname}_{method}_{power}')
    script: 'code/plot_cost_map.py'


rule plot_expenditure_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/expenditure_maps_{nname}_{method}_{power}')
    script: 'code/plot_expenditure_map.py'



rule plot_bars:
    input:
        network = 'resources/{nname}.nc',
        payments = 'resources/payments_{nname}_{method}_{power}.nc'
    output: 'figures/bars_{nname}_{method}_{power}.png'
    script: 'code/plot_bars.py'


# Local Variables:
# mode: python
# End:
