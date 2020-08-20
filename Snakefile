configfile: "config.yaml"


wildcard_constraints:
    clusters="[0-9][0-9]",
    field="(bf|gf)"

rule all:
    input:
        expand('figures/{nname}/network.png', **config['analysis']),
        expand('figures/{nname}/total_costs.png', **config['analysis']),
        expand('figures/{nname}/average_price.png', **config['analysis']),
        expand('figures/{nname}/average_demand.png', **config['analysis']),
        expand('figures/{nname}/power_mix_{method}_{power}.png', **config['analysis']),
        expand('figures/{nname}/bars_{method}_{power}.png', **config['analysis']),
        expand('figures/{nname}/locality_{method}_{power}.png', **config['analysis']),
        expand('figures/{nname}/maps_expenditure_{method}_{power}', **config['analysis']),
        expand('figures/{nname}/maps_price_{method}_{power}', **config['analysis']),
        expand('figures/{nname}/maps_sparcity_expenditure_{method}_{power}', **config['analysis']),


rule test:
    input:
        expand('figures/{nname}/network.png', **config['test']),
        expand('figures/{nname}/total_costs.png', **config['test']),
        expand('figures/{nname}/average_price.png', **config['test']),
        expand('figures/{nname}/average_demand.png', **config['test']),
        expand('figures/{nname}/power_mix_{method}_{power}.png', **config['test']),
        expand('figures/{nname}/bars_{method}_{power}.png', **config['test']),
        expand('figures/{nname}/locality_{method}_{power}.png', **config['test']),
        expand('figures/{nname}/maps_expenditure_{method}_{power}', **config['test']),
        expand('figures/{nname}/maps_price_{method}_{power}', **config['test']),
        expand('figures/{nname}/maps_sparcity_expenditure_{method}_{power}', **config['test']),


subworkflow pypsade:
    workdir: "pypsa-de"
    snakefile: "pypsa-de/Snakefile"
    configfile: "config_germany.yaml"

subworkflow pypsaeur:
    workdir: "pypsa-eur"
    snakefile: "pypsa-eur/Snakefile"
    configfile: "config_europe.yaml"


rule solve_and_sanitize_german_network:
    input:
        network = pypsade('networks/elec_s_{clusters}_ec_lvopt_Ep.nc'),
        regions = pypsade('resources/regions_onshore_elec_s_{clusters}.geojson')
    output: 
        network = 'resources/de{clusters}{field}.nc',
        regions = 'resources/de{clusters}{field}_regions.geojson'
    script: 'code/solve_and_sanitize_network.py'


rule solve_and_sanitize_european_network:
    input:
        network = pypsaeur('networks/elec_s_{clusters}_ec_lvopt_Ep.nc'),
        regions = pypsaeur('resources/regions_onshore_elec_s_{clusters}.geojson')
    output: 
        network = 'resources/eur{clusters}{field}.nc',
        regions = 'resources/eur{clusters}{field}_regions.geojson'
    script: 'code/solve_and_sanitize_network.py'


rule solve_tiny_network:
    output: 'resources/acdc.nc'
    script: 'code/solve_tiny_network.py'


rule create_test_network:
    input: 
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson'
    output: 
        network = 'resources/test-{nname}.nc',
        regions = 'resources/test-{nname}_regions.geojson'
    script: 'code/create_test_network.py'

rule check_shadowprices:
    input:
        network = 'resources/{nname}.nc'
    output: 
        touch('resources/{nname}_shadowprices_checked')
    script: 'code/check_shadowprices.py'


rule allocate_network:
    input:
        network = 'resources/{nname}.nc',
        check = 'resources/{nname}_shadowprices_checked'
    output:
        costs = 'resources/costs_{nname}_{method}_{power}.nc',
        payments = 'resources/payments_{nname}_{method}_{power}.nc', 
        sparcity = 'resources/sparcity_costs_{nname}_{method}_{power}.nc',
        power_mix = 'resources/power_mix_{nname}_{method}_{power}.nc'
    script: 'code/allocate_network.py'


rule plot_network:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson'
    output: 'figures/{nname}/network.png'
    script: 'code/plot_network.py'


rule plot_total_costs:
    input:
        costs = 'resources/costs_{nname}_ptpf_net.nc',
        network = 'resources/{nname}.nc',
    output: 'figures/{nname}/total_costs.png'
    script: 'code/plot_total_costs.py'


rule plot_power_mix:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        power_mix = 'resources/power_mix_{nname}_{method}_{power}.nc'
    output: 'figures/{nname}/power_mix_{method}_{power}.png'
    script: 'code/plot_power_mix.py'


rule plot_total_cost_comparison:
    input:
        costs_gf = 'resources/costs_{nname_wo_field}gf_ptpf_net.nc',
        costs_bf = 'resources/costs_{nname_wo_field}bf_ptpf_net.nc',
        network_gf = 'resources/{nname_wo_field}gf.nc',
        network_bf = 'resources/{nname_wo_field}bf.nc'
    output: 'figures/total_costs_{nname_wo_field}.png'
    script: 'code/plot_total_cost_comparison.py'


rule plot_average_price:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
    output: 'figures/{nname}/average_price.png'
    script: 'code/plot_map.py'

rule plot_average_demand:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
    output: 'figures/{nname}/average_demand.png'
    script: 'code/plot_map.py'


rule plot_expenditure_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/{nname}/maps_expenditure_{method}_{power}')
    script: 'code/plot_maps_transfer.py'

rule plot_price_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/{nname}/maps_price_{method}_{power}')
    script: 'code/plot_maps_transfer.py'

rule plot_sparcity_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/sparcity_costs_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/{nname}/maps_sparcity_expenditure_{method}_{power}')
    script: 'code/plot_maps_transfer.py'


rule plot_locality:
    input:
        network = 'resources/{nname}.nc',
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output: 'figures/{nname}/locality_{method}_{power}.png'
    script: 'code/plot_locality.py'


rule plot_bars:
    input:
        network = 'resources/{nname}.nc',
        payments = 'resources/payments_{nname}_{method}_{power}.nc'
    output: 'figures/{nname}/bars_{method}_{power}.png'
    script: 'code/plot_bars.py'


rule plot_allocated_payment:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output:
        'figures/{nname}/allocated_payment_sink_{sink}.png'



# Local Variables:
# mode: python
# End:
