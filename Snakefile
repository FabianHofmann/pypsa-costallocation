configfile: "config.yaml"


wildcard_constraints:
    clusters="\d+",
    field="(bf|gf)",
    method="(ptpf|ebe)",
    expenditure="(all|capex|opex)",


rule all:
    input:
        expand('figures/{nname}/network.png', **config['analysis']),
        expand('figures/{nname}/total_costs.png', **config['analysis']),
        expand('figures/{nname}/average_price.png', **config['analysis']),
        expand('figures/{nname}/average_demand.png', **config['analysis']),
        expand('figures/{nname}/subsidy.png', **config['analysis']),
        expand('figures/{nname}/capex_duration_curve.png', **config['analysis']),
        expand('figures/{nname}/operation_high_expenditure.png', **config['analysis']),
        expand('figures/{nname}/power_mix_{method}_{power}.png', **config['analysis']),
        expand('figures/{nname}/bars_{method}_{power}.png', **config['analysis']),
        expand('figures/{nname}/locality_{expenditure}_{method}_{power}.png', **config['analysis']),
        expand('figures/{nname}/maps_expenditure_{method}_{power}', **config['analysis']),
        expand('figures/{nname}/maps_price_{method}_{power}', **config['analysis']),
        expand('figures/{nname}/maps_scarcity_expenditure_{method}_{power}', **config['analysis']),
        expand('figures/{nname}/maps_scarcity_price_{method}_{power}', **config['analysis']),
        expand('figures/{nname}/{method}_{power}_to_{sink}.png', **config['analysis']),
        expand('tables/{nname}_{method}_{power}', **config['analysis']),



rule test:
    input:
        expand('figures/{nname}/network.png', **config['test']),
        expand('figures/{nname}/total_costs.png', **config['test']),
        expand('figures/{nname}/average_price.png', **config['test']),
        expand('figures/{nname}/average_demand.png', **config['test']),
        expand('figures/{nname}/subsidy.png', **config['test']),
        expand('figures/{nname}/capex_duration_curve.png', **config['test']),
        expand('figures/{nname}/operation_high_expenditure.png', **config['test']),
        expand('figures/{nname}/power_mix_{method}_{power}.png', **config['test']),
        expand('figures/{nname}/bars_{method}_{power}.png', **config['test']),
        expand('figures/{nname}/locality_{expenditure}_{method}_{power}.png', **config['test']),
        expand('figures/{nname}/maps_expenditure_{method}_{power}', **config['test']),
        expand('figures/{nname}/maps_price_{method}_{power}', **config['test']),
        expand('figures/{nname}/maps_scarcity_expenditure_{method}_{power}', **config['test']),
        expand('figures/{nname}/maps_scarcity_price_{method}_{power}', **config['test']),
        expand('figures/{nname}/{method}_{power}_to_{sink}.png', **config['test']),
        expand('tables/{nname}_{method}_{power}', **config['test']),


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
        network = 'resources/de{clusters}{field}.nc',
        regions = 'resources/de{clusters}{field}_regions.geojson'
    script: 'code/solve_and_sanitize_network.py'


rule solve_and_sanitize_european_network:
    input:
        network = pypsaeur('networks/elec_s_{clusters}_ec_lv1.2_Ep.nc'),
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
        scarcity = 'resources/scarcity_costs_{nname}_{method}_{power}.nc',
        power_mix = 'resources/power_mix_{nname}_{method}_{power}.nc'
    script: 'code/allocate_network.py'

rule plot_graphical_abstract:
    input:
        network = 'resources/de50bf.nc',
        regions = 'resources/de50bf_regions.geojson',
    output: 'figures/graphical_abstract.png'
    script: 'code/graphical_abstract.py'


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

rule plot_scarcity_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/scarcity_costs_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/{nname}/maps_scarcity_expenditure_{method}_{power}')
    script: 'code/plot_maps_transfer.py'

rule plot_scarcity_price_maps:
    input:
        network = 'resources/{nname}.nc',
        regions = 'resources/{nname}_regions.geojson',
        costs = 'resources/scarcity_costs_{nname}_{method}_{power}.nc'
    output:
        folder = directory('figures/{nname}/maps_scarcity_price_{method}_{power}')
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
        costs = 'resources/costs_{nname}_{method}_{power}.nc'
    output: 'figures/{nname}/locality_{expenditure}_{method}_{power}.png'
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
        'figures/{nname}/{method}_{power}_to_{sink}.png'
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
        costs = 'resources/costs_{nname}_{method}_{power}.nc',
        scarcity = 'resources/scarcity_costs_{nname}_{method}_{power}.nc',
    output:
        folder = directory('tables/{nname}_{method}_{power}')
    script: 'code/write_out_tables.py'


# Local Variables:
# mode: python
# End:
