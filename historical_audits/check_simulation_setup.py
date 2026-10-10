#!/usr/bin/env python3
"""Inspect one existing temperature-quench setup without submitting or running MD.

Standard-library only. Run in the same environment intended for the pilot.
This checks prerequisites; it cannot validate an MD protocol scientifically.
"""
import argparse
import importlib
import importlib.util
import json
import shutil
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--group', choices=['YU_E', 'NU_E', 'NU_NE'], default='YU_E')
    parser.add_argument('--protein', default='P16152')
    args = parser.parse_args()
    folder = args.repo / '1_Simulations_and_analyses' / '1_Temperature_quench' / args.group / args.protein
    path = folder / 'Tq.cntrl'
    config = {}
    for raw in path.read_text().splitlines():
        line = raw.split('#', 1)[0].strip()
        if '=' in line:
            key, value = line.split('=', 1)
            config[key.strip()] = value.strip()
    files = {}
    for key in ['psf', 'top', 'param', 'starting_strucs', 'secondary_structure_def']:
        value = config.get(key)
        files[key] = {'configured_path': value, 'exists': bool(value) and (folder / value).is_file()}
    dependencies = {}
    for name in ['numpy', 'pandas', 'mdtraj', 'parmed', 'pyemma', 'openmm', 'simtk']:
        dependencies[name] = importlib.util.find_spec(name) is not None
    engines = []
    engine_error = None
    for name in ['openmm', 'simtk.openmm']:
        try:
            mod = importlib.import_module(name)
            engines = [mod.Platform.getPlatform(i).getName() for i in range(mod.Platform.getNumPlatforms())]
            break
        except ImportError as error:
            engine_error = str(error)
        except Exception as error:
            engine_error = f'{type(error).__name__}: {error}'
    commands = {name: shutil.which(name) for name in [
        'sbatch', 'opt_nscal.py', 'temperature_quenching.py',
        'calc_native_contact_fraction.pl', 'calc_entanglement_number_v2.pl',
        'calc_chirality_number.pl', 'stride']}
    result = {'selected_protein': args.protein, 'original_group': args.group,
              'python': sys.version.split()[0], 'platform': sys.platform,
              'setup_files': files, 'configuration': config,
              'python_packages_found': dependencies, 'executables': commands,
              'openmm_platforms': engines, 'openmm_import_error': engine_error if not engines else None,
              'raw_trajectory_files_in_selected_folder': len(list(folder.rglob('*.dcd'))),
              'notes': [
                  'No notebook, simulation, analysis wrapper, or SLURM submission was run.',
                  'Existing labels concern ubiquitination; assign ED/NED separately for the new project.',
                  'Supplied defaults are production-scale. A fresh pilot needs a copied directory and reviewed run length, trajectory count, restart behavior, and job resources.',
                  'Analysis wrappers/parsers also hard-code trajectory counts. Update all relevant copies together.',
                  'Finding files and executables does not establish force-field validity or successful MD integration.']}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
