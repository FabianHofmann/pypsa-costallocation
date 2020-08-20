configfile: "config.yaml"


wildcard_constraints:
    clusters="[0-9][0-9]",
    field="(bf|gf)"

rule all:
    input:
        expand('figures/{nname}.png', **config['analysis']),
        expand('figures/total_costs_{nname}.png', **config['analysis']),
        expand('figures/average_price_{nname}.png', **config['analysis']),
        expand('figures/average_demand_{nname}.png', **config['analysis']),
        expand('figures/power_mix_{nname}_{method}_{power}.png', **config['analysis']),
        expand('figures/bars_{nname}_{method}_{power}.png', **config['analysis']),
        expand('figures/maps_expenditure_{nname}_{method}_{power}', **config['analysis']),
        expand('figures/maps_price_{nname}_{method}_{power}', **config['analysis'])


rule test:
    input:
        expand('figures/{nname}.png', **config['test']),
        expand('figures/total_costs_{nname}.png', **config['test']),
        expand('figures/average_price_{nname}.png', **config['test']),
        expand('figures/average_demand_{nname}.png', **config['test']),
        expand('figures/power_mix_{nname}_{method}_{power}.png', **config['test']),
        expand('figures/bars_{nname}_{method}_{power}.png', **config['test']),
        expand('figures/maps_expenditure_{nname}_{method}_{power}', **config['test']),
        expand('figures/maps_price_{nname}_{method}_{power}', **config['test'])


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
        sparcity = 'resources/sparcity_cost_{nname}_{method}_{power}.nc',
        power_mix = 'resources/power_mix_{nname}_{method}_{power}.nc'
    script: 'code/allocate_network.py'


rule plot_network:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson'
    output: 'figures/{nname}.png'
    script: 'code/plot_network.py'


rule plot_total_costs:
    input:
        costs = 'resources/costs_{nname}_ptpf_net.nc',
        network = 'resources/{nname}.nc',
    output: 'figures/total_costs_{nname}.png'
    script: 'code/plot_total_costs.py'


rule plot_power_mix:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        power_mix = 'resources/power_mix_{nname}_{method}_{power}.nc'
    output: 'figures/power_mix_{nname}_{method}_{power}.png'
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
    output: 'figures/average_price_{nname}.png'
    script: 'code/plot_map.py'

rule plot_average_demand:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
    output: 'figures/average_demand_{nname}.png'
    script: 'code/plot_map.py'


# rule plot_payment_maps:
#     input:
#         network = 'resources/{nname}.nc',
#         regions = 'resources/{nname}_regions.geojson',
#         payments = 'resources/payments_{nname}_{method}_{power}.nc'
#     output:
#         folder = directory('figures/maps_payment_{nname}_{method}_{power}')
#     script: 'code/plot_maps_payment.py'


rule plot_expenditure_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/maps_expenditure_{nname}_{method}_{power}')
    script: 'code/plot_maps_expenditure.py'


rule plot_bars:
    input:
        network = 'resources/{nname}.nc',
        payments = 'resources/payments_{nname}_{method}_{power}.nc'
    output: 'figures/bars_{nname}_{method}_{power}.png'
    script: 'code/plot_bars.py'

rule plot_allocated_payment:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output:
        'figures/allocated_payment_sink_{sink}.png'



# Local Variables:
# mode: python
# End:
