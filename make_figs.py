import subprocess

files = [
    'fig_step',
    'fig_spec',
    'fig_delay',
    'fig_cg2_compare',
    'fig_restart_compare',
    'fig_2norm_opt',
    'timing',
]

for file in files:
    subprocess.run(f'jupyter nbconvert --execute --to notebook --inplace --allow-errors --ExecutePreprocessor.timeout=-1 {file}.ipynb',shell=True)

