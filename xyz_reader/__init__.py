"""
Extended XYZ reader
git clone --recursive https://gitlab.com/ase/ase.git

cd ase
python3 -m ase test fio/extxyz.py
"""
from abc import ABCMeta, abstractmethod
from pathlib import Path
from typing import Union

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from xyz_reader.tools import key_val_str_to_dict


class XYZReader(metaclass=ABCMeta):
    def __init__(self, file: Union[str, Path]):
        self.file = file
        self._fp = None

    def __enter__(self):
        self._fp = open(self.file, 'r').__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self._fp.__exit__(exc_type, exc_val, exc_tb)

    def __iter__(self):
        while self._fp is not None:
            yield self.read_frame()

    @abstractmethod
    def read_frame(self):
        pass

    @abstractmethod
    def read_atoms(self):
        pass


class SimpleXYZ(XYZReader):

    def read_frame(self):
        n_atoms = int(next(self._fp))
        comments = next(self._fp).rstrip()
        species, positions = self._read_atoms_properties(n_atoms)

        # TODO: VEC ....

        return n_atoms, comments, species, positions

    def _read_atoms_properties(self, n_atoms: int):

        species, positions = [], []

        for _ in range(n_atoms):
            data = next(self._fp).split()

            species.append(int(data[0]))
            positions.append([float(data[1]), float(data[2]), float(data[3])])

        return species, positions

    def read_atoms(self):
        for n_atoms, comments, species, positions in self:
            atoms = Atoms(numbers=species, positions=positions)
            atoms.info.update({'comment': comments})

            yield atoms


class ExtendedXYZ(XYZReader):
    DTYPES = {
        'R': float,
        'I': int,
        'S': str,
        'L': lambda x: True if x in ['T', 'True'] else False
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._cache = {}

    def read_frame(self):
        n_atoms = int(next(self._fp))
        frame_properties, property_types = self._read_frame_properties()
        atoms_properties = self._read_atoms_properties(n_atoms, property_types)

        return n_atoms, frame_properties, atoms_properties

    def _read_frame_properties(self):
        line = next(self._fp)
        # data = {key: single_value or quoted_value for key, single_value, quoted_value in self.REGEXP.findall(line)}
        # frame_properties = {key: convert(value) for key, value in data.items()}

        frame_properties = key_val_str_to_dict(line)
        property_types = frame_properties.pop('Properties', None)

        return frame_properties, property_types

    def _get_converter(self, property_types):
        converter = self._cache.get(property_types, None)
        if converter is None:
            # Format is "[NAME:TYPE:NCOLS]...]", e.g. "species:S:1:pos:R:3".
            props = property_types.split(':')
            converter = tuple((name, self.DTYPES[dtype], int(ncols))
                              for name, dtype, ncols
                              in zip(props[0::3], props[1::3], props[2::3]))
            self._cache.update({property_types: converter})
        return converter

    def _read_atoms_properties(self, n_atoms: int, property_types: str):

        converter = self._get_converter(property_types)
        results = {key: [] for key, _, _ in converter}

        for _ in range(n_atoms):
            data = next(self._fp).split()
            ind = 0
            for (name, dtype, ncols) in converter:
                if ncols == 1:
                    results[name].append(data[ind])
                    ind += ncols
                    continue

                results[name].append([dtype(x) for x in data[ind:ind + ncols]])
                ind += ncols

        return results

    def read_atoms(self, energy_label='energy', forces_label='forces'):
        for n_atoms, frame_properties, atoms_properties in self:

            cell = frame_properties.pop('Lattice', None)
            if cell is not None:
                cell = np.array(cell, order='F').reshape((3, 3))

            pbc = frame_properties.pop('pbc', None)
            if cell is not None and pbc is None:
                pbc = [True, True, True]

            symbols = atoms_properties.pop('species', None)
            numbers = atoms_properties.pop('Z', None)
            if symbols is not None:
                symbols = [s.capitalize() for s in symbols]
                numbers = None

            positions = atoms_properties.pop('pos', None)
            charges = atoms_properties.pop('charge', None)

            atoms = Atoms(symbols=symbols,
                          positions=positions,
                          numbers=numbers,
                          charges=charges,
                          cell=cell,
                          pbc=pbc)

            # Load results of previous calculations into SinglePointCalculator
            results = {}

            energy = frame_properties.pop(energy_label, None)
            if energy is not None:
                results['energy'] = energy

            magmom = frame_properties.pop('magmom', None)
            if magmom is not None:
                results['magmom'] = magmom

            free_energy = frame_properties.pop('free_energy', None)
            if free_energy is not None:
                results['free_energy'] = free_energy

            forces = atoms_properties.pop(forces_label, None)
            if forces is not None:
                results['forces'] = forces

            dipole = atoms_properties.pop('dipole', None)
            # TODO: Make sure that it has the proper representation
            if dipole is not None:
                results['dipole'] = dipole

            charges = atoms_properties.pop('charges', None)
            # TODO: Make sure that it has the proper representation
            if charges is not None:
                results['charges'] = charges

            magmoms = atoms_properties.pop('magmoms', None)
            # TODO: Make sure that it has the proper representation
            if magmoms is not None:
                results['magmoms'] = magmoms

            stress = atoms_properties.pop('stress', None)
            if stress is not None:
                stress = np.array(stress).reshape((3, 3), order='F')
                stress = np.array([stress[0, 0],
                                   stress[1, 1],
                                   stress[2, 2],
                                   stress[1, 2],
                                   stress[0, 2],
                                   stress[0, 1]])
                results['stress'] = stress

            if results:
                calculator = SinglePointCalculator(atoms, **results)
                atoms.set_calculator(calculator)

            # Storing all the remaining properties in the info
            atoms.info.update(frame_properties)
            atoms.info.update(atoms_properties)

            yield atoms


# create aliases for read/write functions
def read_extxyz(*args, **kwargs):
    pass


def write_extxyz(*args, **kwargs):
    pass


if __name__ == '__main__':
    from pathlib import Path
    import json

    file = Path('../test/complex.xyz')

    # file = directory / 'GAP_6.xyz'

    with ExtendedXYZ(file) as reader:
        for atoms in reader.read_atoms(forces_label='force'):
            print('==========================')
            print(atoms)

    with ExtendedXYZ(file) as reader:
        for n_atoms, frame_properties, atoms_properties in reader:
            # print('==========================')
            print(f'natoms: {n_atoms}')
            print(frame_properties)
            # print(json.dumps(frame_properties, indent=2))
            print(atoms_properties)
