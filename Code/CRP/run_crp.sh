#!/bin/bash
#make sure to make the output subject folder
echo running on Spat31
python crp_pipeline_runner.py -i 'input_files/Spat31/inp.csv' -c '7'
echo Done with Spat 31
echo Running on Epat31
python crp_pipeline_runner.py -i 'input_files/Epat31/inp.csv' -c '7'
echo Done with Epat31
echo Running on Epat27
python crp_pipeline_runner.py -i 'input_files/Epat27/inp.csv' -c '7'
echo Done with Epat27
echo Running on Epa26
python crp_pipeline_runner.py -i 'input_files/Epat26/inp.csv' -c '7'
echo Done with Epat26!
echo running on Epat30
python crp_pipeline_runner.py -i 'input_files/Epat30/inp.csv' -c '7'
echo Done with Epat30
echo running on Epat35
python crp_pipeline_runner.py -i 'input_files/Epat35/inp.csv' -c '7'
echo done with Epat25
echo running on Epat37
python crp_pipeline_runner.py -i 'input_files/Epat37/inp.csv' -c '7'
echo Done running Epat37
echo running on Epat38
python crp_pipeline_runner.py -i 'input_files/Epat38/inp.csv' -c '7'
echo Done running on Epat38
echo Running on Epat39
python crp_pipeline_runner.py -i 'input_files/Epat39/inp.csv' -c '7'
echo Done running Epat 39
echo running Epat43
python crp_pipeline_runner.py -i 'input_files/Epat43/inp.csv' -c '7'
echo Done with Epat43
python crp_pipeline_runner.py -i 'input_files/Spat30/inp.csv' -c '7'
echo Done with Spat30

python crp_pipeline_runner.py -i 'input_files/Spat34/inp.csv' -c '7'
echo done with Spat34
python crp_pipeline_runner.py -i 'input_files/Spat36/inp.csv' -c '7'
echo done with Spat36
python crp_pipeline_runner.py -i 'input_files/Spat37/inp.csv' -c '7'
echo done with Spat37
python crp_pipeline_runner.py -i 'input_files/Spat52/inp.csv' -c '7'
echo done with Spat52