from yutility import orbitals, squeeze_list
import unittest
import os



## first test

class ClosedShellSinglePoint(unittest.TestCase):
    def setUp(self):
        self.orbs = orbitals.Orbitals('rkf/substrate_cat_complex.rkf', moleculename='Substrate')
        self.lumo = self.orbs.mos['LUMO']
        self.homo = self.orbs.mos['HOMO-5']
        self.sfo = squeeze_list(self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1))


    def test_sfo_select(self):
        self.assertEqual(self.orbs.sfos['12A'], squeeze_list(self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1)))

    def test_sfo_select2(self):
        self.assertEqual(self.orbs.sfos.get_sfo(index=12), self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1))

    def test_sfo_select3(self):
        self.assertEqual(self.orbs.sfos['C:1(1P:z)'], squeeze_list(self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1)))

    def test_sfo_select4(self):
        self.assertEqual([self.orbs.sfos['C:1(1P:x)'], self.orbs.sfos['C:1(1P:y)'], self.orbs.sfos['C:1(1P:z)']], self.orbs.sfos.get_sfo(orbname='1P', fragidx=1))


    def test_mo_select(self):
        self.assertEqual(len(self.orbs.mos['HOMO-5':'HOMO']), 6)

    def test_mo_select2(self):
        self.assertEqual(self.orbs.mos['HOMO-5'], self.orbs.mos['77A'])

    def test_mo_select3(self):
        self.assertEqual(self.orbs.mos['HOMO-5'], self.orbs.mos[77])


    def test_lumo_energy(self):
        self.assertTrue(abs(-4.730192 - self.lumo.energy) < 0.0001)

    def test_lumo_name(self):
        self.assertEqual(self.lumo.relative_name, 'LUMO')

    def test_lumo_occupied(self):
        self.assertFalse(self.lumo.occupied)

    def test_lumo_occupation(self):
        self.assertEqual(self.lumo.occupation, 0)

    def test_lumo_offset(self):
        self.assertEqual(self.lumo.relindex, 1)

    def test_lumo_spin(self):
        self.assertEqual(self.lumo.spin, 'AB')

    def test_lumo_symmetry(self):
        self.assertEqual(self.lumo.symmetry, 'A')

    def test_lumo_fragment(self):
        self.assertEqual(self.lumo.moleculename, 'Substrate')

    def test_lumo_rkf_path(self):
        self.assertEqual(self.lumo.kfpath, os.path.abspath('rkf/substrate_cat_complex.rkf'))


    def test_homo_energy(self):
        self.assertTrue(abs(-8.579126358615937 - self.homo.energy) < 0.0001)

    def test_homo_name(self):
        self.assertEqual(self.homo.relative_name, 'HOMO-5')

    def test_homo_occupied(self):
        self.assertTrue(self.homo.occupied)

    def test_homo_occupation(self):
        self.assertEqual(self.homo.occupation, 2)

    def test_homo_offset(self):
        self.assertEqual(self.homo.relindex, -5)

    def test_homo_spin(self):
        self.assertEqual(self.homo.spin, 'AB')

    def test_homo_symmetry(self):
        self.assertEqual(self.homo.symmetry, 'A')

    def test_homo_fragment(self):
        self.assertEqual(self.homo.moleculename, 'Substrate')

    def test_homo_rkf_path(self):
        self.assertEqual(self.homo.kfpath, os.path.abspath('rkf/substrate_cat_complex.rkf'))


    def test_sfo_name(self):
        self.assertEqual(self.sfo.name, '1P:z')

    def test_sfo_occupation(self):
        self.assertEqual(round(self.sfo.occupation, 2), 0.67)

    def test_sfo_spin(self):
        self.assertEqual(self.sfo.spin, 'AB')

    def test_sfo_symmetry(self):
        self.assertEqual(self.sfo.symmetry, 'A')

    def test_sfo_fragment(self):
        self.assertEqual(self.sfo.fragment, 'C')

    def test_sfo_fragmentindex(self):
        self.assertEqual(self.sfo.fragment_index, 1)

    def test_sfo_rkf_path(self):
        self.assertEqual(self.sfo.kfpath, os.path.abspath('rkf/substrate_cat_complex.rkf'))


    def test_coeff(self):
        self.assertEqual(round(self.sfo @ self.lumo, 4), -0.4617)

    def test_coeff2(self):
        self.assertEqual(round(self.lumo @ self.sfo, 4), -0.4617)

    def test_overlap(self):
        self.assertEqual(round(self.orbs.sfos['268A'] @ self.orbs.sfos['12A'], 4), -0.0092)

    def test_overlap2(self):
        self.assertEqual(round(self.orbs.sfos['12A'] @ self.orbs.sfos['268A'], 4), -0.0092)



class ClosedShellEDA(unittest.TestCase):
    def setUp(self):
        self.orbs = orbitals.Orbitals('rkf/BH3NH3.rkf', moleculename='LAadduct')
        self.lumo = self.orbs.mos['LUMO+2']
        self.homo = self.orbs.mos['HOMO']
        self.sfo = self.orbs.sfos['Donor(HOMO-3)']

    def test_unrestricted_sfo(self):
        self.assertFalse(self.orbs.sfos.is_unrestricted)

    def test_unrestricted_mo(self):
        self.assertFalse(self.orbs.mos.is_unrestricted)

    def test_relativistic_sfo(self):
        self.assertTrue(self.orbs.sfos.is_relativistic)

    def test_relativistic_mo(self):
        self.assertTrue(self.orbs.mos.is_relativistic)

    def test_atomic_sfo(self):
        self.assertFalse(self.orbs.sfos.uses_atomic_fragments)


    def test_sfo_select(self):
        self.assertEqual(set(self.orbs.sfos['12A']), {self.orbs.sfos['Acceptor(12A)'], self.orbs.sfos['Donor(12A)']})

    def test_sfo_select2(self):
        self.assertEqual(squeeze_list(self.orbs.sfos.get_sfo(index=12)), self.orbs.sfos[12])

    def test_sfo_select3(self):
        self.assertEqual(self.orbs.sfos['Acceptor(5A)'], self.orbs.sfos['Acceptor(LUMO)'])

    def test_sfo_select4(self):
        self.assertEqual(set(self.orbs.sfos['HOMO']), {self.orbs.sfos['Acceptor(HOMO)'], self.orbs.sfos['Donor(HOMO)']})


    def test_mo_select(self):
        self.assertEqual(len(self.orbs.mos['HOMO-5':'HOMO']), 6)

    def test_mo_select2(self):
        self.assertEqual(self.orbs.mos['HOMO-5'], self.orbs.mos['4A'])

    def test_mo_select3(self):
        self.assertEqual(self.orbs.mos['HOMO-5'], self.orbs.mos[4])


    def test_sfo_fragments(self):
        self.assertEqual(self.orbs.sfos.fragments, {'Donor', 'Acceptor'})


    def test_lumo_energy(self):
        self.assertTrue(abs(0.5228762356838158 - self.lumo.energy) < 0.0001)

    def test_lumo_name(self):
        self.assertEqual(self.lumo.relative_name, 'LUMO+2')

    def test_lumo_occupied(self):
        self.assertFalse(self.lumo.occupied)

    def test_lumo_occupation(self):
        self.assertEqual(self.lumo.occupation, 0)

    def test_lumo_offset(self):
        self.assertEqual(self.lumo.relindex, 3)

    def test_lumo_spin(self):
        self.assertEqual(self.lumo.spin, 'AB')

    def test_lumo_symmetry(self):
        self.assertEqual(self.lumo.symmetry, 'A')

    def test_lumo_fragment(self):
        self.assertEqual(self.lumo.moleculename, 'LAadduct')

    def test_lumo_rkf_path(self):
        self.assertEqual(self.lumo.kfpath, os.path.abspath('rkf/BH3NH3.rkf'))


    def test_homo_energy(self):
        self.assertTrue(abs(-6.447661 - self.homo.energy) < 0.0001)

    def test_homo_name(self):
        self.assertEqual(self.homo.relative_name, 'HOMO')

    def test_homo_occupied(self):
        self.assertTrue(self.homo.occupied)

    def test_homo_occupation(self):
        self.assertEqual(self.homo.occupation, 2)

    def test_homo_offset(self):
        self.assertEqual(self.homo.relindex, 0)

    def test_homo_spin(self):
        self.assertEqual(self.homo.spin, 'AB')

    def test_homo_symmetry(self):
        self.assertEqual(self.homo.symmetry, 'A')

    def test_homo_fragment(self):
        self.assertEqual(self.homo.moleculename, 'LAadduct')

    def test_homo_rkf_path(self):
        self.assertEqual(self.homo.kfpath, os.path.abspath('rkf/BH3NH3.rkf'))


    def test_sfo_name(self):
        self.assertEqual(self.sfo.name, '2A')

    def test_sfo_occupation(self):
        self.assertEqual(self.sfo.occupation, 2)

    def test_sfo_indexname(self):
        self.assertEqual(self.sfo.index_name, '2A')

    def test_sfo_relname(self):
        self.assertEqual(self.sfo.relative_name, 'Donor(HOMO-3)')

    def test_sfo_spin(self):
        self.assertEqual(self.sfo.spin, 'AB')

    def test_sfo_symmetry(self):
        self.assertEqual(self.sfo.symmetry, 'A')

    def test_sfo_fragment(self):
        self.assertEqual(self.sfo.fragment, 'Donor')

    def test_sfo_fragmentindex(self):
        self.assertEqual(self.sfo.fragment_index, 1)

    def test_sfo_rkf_path(self):
        self.assertEqual(self.sfo.kfpath, os.path.abspath('rkf/BH3NH3.rkf'))


    def test_coeff(self):
        sfo = self.orbs.sfos['Acceptor(HOMO-2)']
        mo = self.orbs.mos['HOMO-2']
        self.assertEqual(round(sfo @ mo, 4), -0.7937)

    def test_coeff2(self):
        sfo = self.orbs.sfos['Donor(LUMO+1)']
        mo = self.orbs.mos['LUMO+1']
        self.assertEqual(round(sfo @ mo, 4), -0.9474)

    def test_overlap(self):
        sfo1 = self.orbs.sfos['Acceptor(2A)']
        sfo2 = self.orbs.sfos['Donor(5A)']
        self.assertEqual(round(sfo1 @ sfo2, 4), -0.3262)

    def test_overlap2(self):
        sfo1 = self.orbs.sfos['Acceptor(LUMO)']
        sfo2 = self.orbs.sfos['Donor(HOMO)']
        self.assertEqual(round(sfo1 @ sfo2, 4), 0.3530)



class OpenShellEDA(unittest.TestCase):
    def setUp(self):
        self.orbs = orbitals.Orbitals('rkf/ethaneEDA.rkf', moleculename='Ethane')
        self.orbs.rename_fragments(('Region_1', 'Region_2'), ('Left', 'Right'))
        self.lumo = self.orbs.mos['LUMO+2_A']
        self.homo = self.orbs.mos['HOMO_B']
        self.sfo = self.orbs.sfos['Left(HOMO-3)_A']

    def test_unrestricted_sfo(self):
        self.assertTrue(self.orbs.sfos.is_unrestricted)

    def test_unrestricted_mo(self):
        self.assertTrue(self.orbs.mos.is_unrestricted)

    def test_relativistic_sfo(self):
        self.assertFalse(self.orbs.sfos.is_relativistic)

    def test_relativistic_mo(self):
        self.assertFalse(self.orbs.mos.is_relativistic)

    def test_atomic_sfo(self):
        self.assertFalse(self.orbs.sfos.uses_atomic_fragments)


    def test_sfo_select(self):
        self.assertEqual(set(self.orbs.sfos['5A']), set(self.orbs.sfos['5A_A', '5A_B']))

    def test_sfo_select2(self):
        self.assertEqual(set(self.orbs.sfos['HOMO-3']), set(self.orbs.sfos['HOMO-3_A', 'HOMO-3_B']))

    def test_sfo_select3(self):
        self.assertEqual(set(self.orbs.sfos['Left(HOMO-3)']), set(self.orbs.sfos['Left(HOMO-3)_A', 'Left(HOMO-3)_B']))

    def test_sfo_select4(self):
        self.assertEqual(self.orbs.sfos.get_sfo(index=12), self.orbs.sfos[12])

    def test_sfo_select5(self):      
        self.assertEqual(set(self.orbs.sfos['HOMO']), set(self.orbs.sfos['Left(HOMO)_A', 'Left(HOMO)_B', 'Right(HOMO)_A', 'Right(HOMO)_B']))

    def test_sfo_select6(self):      
        self.assertEqual(set(self.orbs.sfos['Left(LUMO+19']), set(self.orbs.sfos['Left(LUMO+19)_A', 'Left(LUMO+19)_B']))


    def test_mo_select(self):
        self.assertEqual(len(self.orbs.mos['HOMO-5':'HOMO']), 12)

    def test_mo_select2(self):
        self.assertEqual(self.orbs.mos['HOMO'], self.orbs.mos['9A'])

    def test_mo_select3(self):
        self.assertEqual(self.orbs.mos['HOMO-5'], self.orbs.mos[4])


    def test_sfo_fragments(self):
        self.assertEqual(self.orbs.sfos.fragments, {'Left', 'Right'})


    def test_lumo_energy(self):
        self.assertTrue(abs(0.873039 - self.lumo.energy) < 0.0001)

    def test_lumo_name(self):
        self.assertEqual(self.lumo.name, '12A')

    def test_lumo_full_name(self):
        self.assertEqual(self.lumo.full_name, '12A_A')

    def test_lumo_occupied(self):
        self.assertFalse(self.lumo.occupied)

    def test_lumo_occupation(self):
        self.assertEqual(self.lumo.occupation, 0)

    def test_lumo_relative_name(self):
        self.assertEqual(self.lumo.relative_name, 'LUMO+2_A')

    def test_lumo_offset(self):
        self.assertEqual(self.lumo.relindex, 3)

    def test_lumo_spin(self):
        self.assertEqual(self.lumo.spin, 'A')

    def test_lumo_symmetry(self):
        self.assertEqual(self.lumo.symmetry, 'A')

    def test_lumo_fragment(self):
        self.assertEqual(self.lumo.moleculename, 'Ethane')

    def test_lumo_rkf_path(self):
        self.assertEqual(self.lumo.kfpath, os.path.abspath('rkf/ethaneEDA.rkf'))


    def test_homo_energy(self):
        self.assertTrue(abs(-8.124029 - self.homo.energy) < 0.0001)

    def test_homo_name(self):
        self.assertEqual(self.homo.name, '9A')

    def test_homo_full_name(self):
        self.assertEqual(self.homo.full_name, '9A_B')

    def test_homo_occupied(self):
        self.assertTrue(self.homo.occupied)

    def test_homo_occupation(self):
        self.assertEqual(self.homo.occupation, 1)

    def test_homo_offset(self):
        self.assertEqual(self.homo.relindex, 0)

    def test_homo_spin(self):
        self.assertEqual(self.homo.spin, 'B')

    def test_homo_symmetry(self):
        self.assertEqual(self.homo.symmetry, 'A')

    def test_homo_fragment(self):
        self.assertEqual(self.homo.moleculename, 'Ethane')

    def test_homo_rkf_path(self):
        self.assertEqual(self.homo.kfpath, os.path.abspath('rkf/ethaneEDA.rkf'))


    def test_sfo_name(self):
        self.assertEqual(self.sfo.name, '2A')

    def test_sfo_full_name(self):
        self.assertEqual(self.sfo.full_name, 'Left(2A)_A')

    def test_sfo_relative_name(self):
        self.assertEqual(self.sfo.relative_name, 'Left(HOMO-3)_A')

    def test_sfo_occupation(self):
        self.assertEqual(self.sfo.occupation, 1)

    def test_sfo_spin(self):
        self.assertEqual(self.sfo.spin, 'A')

    def test_sfo_symmetry(self):
        self.assertEqual(self.sfo.symmetry, 'A')

    def test_sfo_fragment(self):
        self.assertEqual(self.sfo.fragment, 'Left')

    def test_sfo_fragmentindex(self):
        self.assertEqual(self.sfo.fragment_index, 2)

    def test_sfo_rkf_path(self):
        self.assertEqual(self.sfo.kfpath, os.path.abspath('rkf/ethaneEDA.rkf'))

    def test_sfo_relindex(self):
        self.assertEqual(self.sfo.relindex, -3)


    # def test_coeff(self):
    #     sfo = self.orbs.sfos['Acceptor(HOMO-2)']
    #     mo = self.orbs.mos['HOMO-2']
    #     self.assertEqual(round(sfo @ mo, 4), -0.7937)

    # def test_coeff2(self):
    #     sfo = self.orbs.sfos['Donor(LUMO+1)']
    #     mo = self.orbs.mos['LUMO+1']
    #     self.assertEqual(round(sfo @ mo, 4), -0.9474)

    # def test_overlap(self):
    #     sfo1 = self.orbs.sfos['Acceptor(2A)']
    #     sfo2 = self.orbs.sfos['Donor(5A)']
    #     self.assertEqual(round(sfo1 @ sfo2, 4), -0.3262)

    # def test_overlap2(self):
    #     sfo1 = self.orbs.sfos['Acceptor(LUMO)']
    #     sfo2 = self.orbs.sfos['Donor(HOMO)']
    #     self.assertEqual(round(sfo1 @ sfo2, 4), 0.3530)



if __name__ == '__main__':
    unittest.main(verbosity=1)
