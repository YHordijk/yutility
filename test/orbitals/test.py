from yutility import orbitals
import unittest



## first test

class ClosedShell(unittest.TestCase):
    def setUp(self):
        self.orbs = orbitals.Orbitals('rkf/substrate_cat_complex.rkf', moleculename='Substrate')
        self.lumo = self.orbs.fmos['LUMO']
        self.homo = self.orbs.fmos['HOMO-5']
        self.sfo = self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1)


    def test_sfo_select(self):
        self.assertEqual(self.orbs.sfos['12A'], self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1))

    def test_sfo_select2(self):
        self.assertEqual(self.orbs.sfos.get_sfo(index=12), self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1))

    def test_sfo_select3(self):
        self.assertEqual(self.orbs.sfos['C:1(1P:z)'], self.orbs.sfos.get_sfo(orbname='1P:z', fragidx=1))

    def test_sfo_select4(self):
        self.assertEqual([self.orbs.sfos['C:1(1P:x)'], self.orbs.sfos['C:1(1P:y)'], self.orbs.sfos['C:1(1P:z)']], self.orbs.sfos.get_sfo(orbname='1P', fragidx=1))


    def test_fmo_select(self):
        self.assertEqual(len(self.orbs.fmos['HOMO-5':'HOMO']), 6)

    def test_fmo_select2(self):
        self.assertEqual(self.orbs.fmos['HOMO-5'], self.orbs.fmos['77A'])

    def test_fmo_select3(self):
        self.assertEqual(self.orbs.fmos['HOMO-5'], self.orbs.fmos[77])


    def test_lumo_energy(self):
        self.assertTrue(abs(-0.173831 - self.lumo.energy) < 0.0001)

    def test_lumo_name(self):
        self.assertEqual(self.lumo.name, 'LUMO')

    def test_lumo_occupied(self):
        self.assertFalse(self.lumo.occupied)

    def test_lumo_occupation(self):
        self.assertEqual(self.lumo.occupation, 0)

    def test_lumo_AMSlevels_name(self):
        self.assertEqual(self.lumo.AMSlevels_name, '83A')

    def test_lumo_offset(self):
        self.assertEqual(self.lumo.offset, 1)

    def test_lumo_spin(self):
        self.assertEqual(self.lumo.spin, 'AB')

    def test_lumo_symmetry(self):
        self.assertEqual(self.lumo.symmetry, 'A')

    def test_lumo_fragment(self):
        self.assertEqual(self.lumo.fragment, 'Substrate')

    def test_lumo_rkf_path(self):
        self.assertEqual(self.lumo.rkf_path, 'rkf/substrate_cat_complex.rkf')


    def test_homo_energy(self):
        self.assertTrue(abs(-0.315277 - self.homo.energy) < 0.0001)

    def test_homo_name(self):
        self.assertEqual(self.homo.name, 'HOMO-5')

    def test_homo_occupied(self):
        self.assertTrue(self.homo.occupied)

    def test_homo_occupation(self):
        self.assertEqual(self.homo.occupation, 2)

    def test_homo_AMSlevels_name(self):
        self.assertEqual(self.homo.AMSlevels_name, '77A')

    def test_homo_offset(self):
        self.assertEqual(self.homo.offset, -5)

    def test_homo_spin(self):
        self.assertEqual(self.homo.spin, 'AB')

    def test_homo_symmetry(self):
        self.assertEqual(self.homo.symmetry, 'A')

    def test_homo_fragment(self):
        self.assertEqual(self.homo.fragment, 'Substrate')

    def test_homo_rkf_path(self):
        self.assertEqual(self.homo.rkf_path, 'rkf/substrate_cat_complex.rkf')


    def test_sfo_name(self):
        self.assertEqual(self.sfo.name, '1P:z')

    def test_sfo_occupation(self):
        self.assertEqual(round(self.sfo.occupation, 2), 0.67)

    def test_sfo_AMSlevels_name(self):
        self.assertEqual(self.sfo.AMSlevels_name, '1P:z')

    def test_sfo_spin(self):
        self.assertEqual(self.sfo.spin, 'AB')

    def test_sfo_symmetry(self):
        self.assertEqual(self.sfo.symmetry, 'A')

    def test_sfo_fragment(self):
        self.assertEqual(self.sfo.fragment, 'C')

    def test_sfo_fragmentindex(self):
        self.assertEqual(self.sfo.fragmentindex, 1)

    def test_sfo_rkf_path(self):
        self.assertEqual(self.sfo.rkf_path, 'rkf/substrate_cat_complex.rkf')


    def test_coeff(self):
        self.assertEqual(round(self.sfo @ self.lumo, 4), -0.4617)

    def test_coeff2(self):
        self.assertEqual(round(self.lumo @ self.sfo, 4), -0.4617)

    def test_overlap(self):
        self.assertEqual(round(self.orbs.sfos['268A'] @ self.orbs.sfos['12A'], 4), -0.0092)

    def test_overlap2(self):
        self.assertEqual(round(self.orbs.sfos['12A'] @ self.orbs.sfos['268A'], 4), -0.0092)



if __name__ == '__main__':
    unittest.main()
