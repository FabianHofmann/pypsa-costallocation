
git clone https://github.com/PyPSA/pypsa-eur.git
git clone https://github.com/PyPSA/pypsa-eur.git ./pypsa-de
git clone https://github.com/PyPSA/technology-data.git

# cp technology-data/outputs/cost_2030.csv pypsa-eur/data/costs.csv
# cp technology-data/outputs/cost_2030.csv pypsa-de/data/costs.csv

git clone https://github.com/PyPSA/pypsa-eur.git ./pypsa-de
conda install environment.yaml
conda activate costallocation

python 
'''
import powerplantmatching as pm
pm.data.OPSD_VRE_country("DE").powerplant.convert_country_to_alpha2().to_csv("pypsa-de/data/custom_powerplants.csv")
'''