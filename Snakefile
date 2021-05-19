configfile: "config.yaml"


wildcard_constraints:
    clusters="\d+",
    expenditure="(all|capex|opex)",


general_figures = [
    'price_decomposition', 'cost_decomposition', #examples
]

network_figures = [
    'network', 'total_costs', 'average_price', 'average_demand', 'subsidy', 'capex_duration_curve',
    'operation_high_expenditure', 'power_mix', 'bars', 'locality_all', 'allocation_to_highest-lmp' ,
]

network_maps = [
    'maps_expenditure', 'maps_price', 'maps_scarcity_price',  'maps_scarcity_expenditure',
]

rule all:
    input:
        expand('figures/{figure}.png', figure=general_figures),
        expand('figures/{nname}/{figure}.png', figure=network_figures, **config['analysis']),
        expand('figures/{nname}/{directory}', directory=network_maps, **config['analysis']),
        expand('tables/{nname}', **config['analysis']),


rule test:
    input:
        expand('figures/{figure}.png', figure=general_figures),
        expand('figures/{nname}/{figure}.png', figure=network_figures, **config['test']),
        expand('tables/{nname}', **config['test']),


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
        network = pypsade('networks/elec_s_{clusters}_ec_lv1.2_Ep.nc'),
        regions = pypsade('resources/regions_onshore_elec_s_{clusters}.geojson')
    output: 
        network = 'resources/de{clusters}.nc',
        regions = 'resources/de{clusters}_regions.geojson'
    script: 'code/solve_and_sanitize_network.py'


rule solve_and_sanitize_european_network:
    input:
        network = pypsaeur('networks/elec_s_{clusters}_ec_lv1.2_Ep.nc'),
        regions = pypsaeur('resources/regions_onshore_elec_s_{clusters}.geojson')
    output: 
        network = 'resources/eur{clusters}.nc',
        regions = 'resources/eur{clusters}_regions.geojson'
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


rule allocate_revenue:
    input:
        network = 'resources/{nname}.nc',
    output:
        revenues = 'resources/revenues_{nname}.nc',
        payments = 'resources/demand_payments_{nname}.nc', 
    script: 'code/allocate_revenue.py'


rule allocate_network:
    input:
        network = 'resources/{nname}.nc',
        check = 'resources/{nname}_shadowprices_checked'
    output:
        costs = 'resources/costs_{nname}.nc',
        payments = 'resources/payments_{nname}.nc', 
        scarcity = 'resources/scarcity_costs_{nname}.nc',
        power_mix = 'resources/power_mix_{nname}.nc'
    script: 'code/allocate_network.py'



rule plot_graphical_abstract:
    input:
        network = 'resources/de50.nc',
        regions = 'resources/de50_regions.geojson',
    output: 'figures/graphical_abstract.png'
    script: 'code/graphical_abstract.py'


rule plot_price_decomposition:
    output: 
        price_decomposition ='figures/price_decomposition.png',
        cost_decomposition = 'figures/cost_decomposition.png'
    script: 'code/price_decomposition.py'



rule plot_network:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson'
    output: 'figures/{nname}/network.png'
    script: 'code/plot_network.py'


rule plot_total_costs:
    input:
        costs = 'resources/costs_{nname}.nc',
        network = 'resources/{nname}.nc',
    output: 'figures/{nname}/total_costs.png'
    script: 'code/plot_total_costs.py'


rule plot_power_mix:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        power_mix = 'resources/power_mix_{nname}.nc'
    output: 'figures/{nname}/power_mix.png'
    script: 'code/plot_power_mix.py'


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
        costs = 'resources/costs_{nname}.nc'
    output:
        folder = directory('figures/{nname}/maps_expenditure')
    script: 'code/plot_maps_transfer.py'

rule plot_price_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}.nc'
    output:
        folder = directory('figures/{nname}/maps_price')
    script: 'code/plot_maps_transfer.py'

rule plot_scarcity_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/scarcity_costs_{nname}.nc'
    output:
        folder = directory('figures/{nname}/maps_scarcity_expenditure')
    script: 'code/plot_maps_transfer.py'

rule plot_scarcity_price_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/scarcity_costs_{nname}.nc'
    output:
        folder = directory('figures/{nname}/maps_scarcity_price')
    script: 'code/plot_maps_transfer.py'

rule plot_subsidy:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
    output:
        'figures/{nname}/subsidy.png'
    script: 'code/plot_subsidy.py'


rule plot_locality:
    input:
        network = 'resources/{nname}.nc',
        costs = 'resources/costs_{nname}.nc'
    output: 'figures/{nname}/locality_{expenditure}.png'
    script: 'code/plot_locality.py'


rule plot_bars:
    input:
        network = 'resources/{nname}.nc',
        payments = 'resources/payments_{nname}.nc'
    output: 'figures/{nname}/bars.png'
    script: 'code/plot_bars.py'


rule plot_allocated_payment:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/costs_{nname}.nc'
    output:
        'figures/{nname}/allocation_to_{sink}.png'
    script: 
        'code/plot_allocated_payment.py'

rule plot_capex_duration_curve:
    input:
        network = 'resources/{nname}.nc',
    output:
        'figures/{nname}/capex_duration_curve.png'
    script: 
        'code/plot_duration_curve.py'


rule plot_operation_high_expenditure:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
    output:
        'figures/{nname}/operation_high_expenditure.png'
    script: 
        'code/plot_operation.py'


rule write_tables:
    input:
        network = 'resources/{nname}.nc',
        costs = 'resources/costs_{nname}.nc',
        scarcity = 'resources/scarcity_costs_{nname}.nc',
    output:
        folder = directory('tables/{nname}')
    script: 'code/write_out_tables.py'


# Local Variables:
# mode: python
# End:
