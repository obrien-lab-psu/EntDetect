from setuptools import setup, find_packages

setup(
    name='EntDetect',
    version='1.1.1',
    description='Entanglement Detection in Protein Structures',
    packages=find_packages(),
    install_requires=[
        'biopython',
        'numpy',
        'scipy',
        'pandas',
        'MDAnalysis',
        'numba',
        'topoly',
        'geom_median',
        'matplotlib',
        'seaborn',
        'scikit-learn',
        'networkx',
        'pyyaml',
        'tqdm',
    ],
    entry_points={
        'console_scripts': [
            'run_entanglement_identification=scripts.run_entanglement_identification:main',
            'run_OP_CGtrajAnal=scripts.run_OP_CGtrajAnal:main',
            'run_OP_AAtrajAnal=scripts.run_OP_AAtrajAnal:main',
            'run_nonnative_entanglement_clustering=scripts.run_nonnative_entanglement_clustering:main',
            'run_change_resolution=scripts.run_change_resolution:main',
            'run_MSM=scripts.run_MSM:main',
            'run_compare_sim2exp=scripts.run_compare_sim2exp:main',
            'run_population_modeling=scripts.run_population_modeling:main',
            'run_montecarlo=scripts.run_montecarlo:main',
        ],
    },
)