import os
from scm import plams
from yutility import log

def get_calc_info(path):
	reader = plams.KFReader(path)
	ret = {}

	# determine if calculation used relativistic corrections
	# if it did, variable 'escale' will be present in 'SFOs'
	# if it didnt, only variable 'energy' will be present
	ret['relativistic'] = ('SFOs', 'escale') in reader

	# determine if SFOs are unrestricted or not
	ret['unrestricted_sfos'] = ('SFOs', 'energy_B') in reader

	# determine if MOs are unrestricted or not
	symlabels = reader.read('Symmetry', 'symlab').strip().split()
	ret['unrestricted_mos'] = (symlabels[0], 'eps_B') in reader

	return ret



for file in os.listdir():
	if file.endswith('.py'):
		continue
	log.log(file)
	log.log(get_calc_info(file))
	log.log()
