git clone hofmann@fias.uni-frankfurt.de/home/vres/data/scm/costallocation.git

git clone https://github.com/PyPSA/pypsa-eur.git
git clone https://github.com/PyPSA/pypsa-eur.git ./pypsa-de
git clone https://github.com/PyPSA/technology-data.git

# cp technology-data/outputs/cost_2030.csv pypsa-eur/data/costs.csv
# cp technology-data/outputs/cost_2030.csv pypsa-de/data/costs.csv

# conda install environment.yaml

# python 
# '''
# import powerplantmatching as pm
# pm.data.OPSD_VRE_country("DE").to_csv("pypsa-de/data/custom_powerplants.csv")
# '''