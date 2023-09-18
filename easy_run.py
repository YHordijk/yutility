from yutility import molecule, settings, run
import shutil


def optimize(xyzfile, outputfile=None, extra_settings=None, keep_calc_dir=False):
	mol = molecule.load(xyzfile)
	character = mol['flags'].get('character', 'minimum')
	charge = mol['flags'].get('charge', 0)
	spinpol = mol['flags'].get('spinpol', 0)
	level = mol['flags'].get('level', 'OLYP/TZ2P/Good')

	mol['molecule'].identifier_data.update(mol['flags'])

	sett = settings.default(level)
	sett = settings.charge(charge, sett)
	sett = settings.spin_polarization(spinpol, sett)
	if character == 'minimum':
		sett = settings.optimization(sett)
	elif character == 'transitionstate':
		sett = settings.transition_state(sett=sett)

	if extra_settings:
		sett.update(extra_settings)

	res = run.run(mol['molecule'], sett, path='__easy_run')
	opt_mol = res.get_main_molecule()
	for opt_at, at in zip(opt_mol, mol['molecule']):
		at.coords = opt_at.coords

	outputfile = outputfile or xyzfile.removesuffix('.xyz') + '_output.xyz'
	molecule.save(mol['molecule'], outputfile, comment=mol['comment'])

	if not keep_calc_dir:
		shutil.rmtree('__easy_run')

	return mol['molecule']


if __name__ == '__main__':
	for i in range(1, 9):
		optimize(f'/Users/yumanhordijk/PhD/yutility/test/xyz/{i}C.xyz', f'/Users/yumanhordijk/PhD/yutility/test/xyz/{i}C.xyz')	
