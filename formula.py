from yutility import ensure_list


def polynomial(coeffs, variables, unity_var='', mult_sign='*', latex_mode=False, coeff_threshold=None):
	'''
	Function to write a polynomial function based on the coefficients and variable names
	'''
	if latex_mode:
		mult_sign = ''
	formula = ''
	for j, (coeff, var) in enumerate(zip(coeffs, variables)):
		if coeff_threshold and coeff < coeff_threshold:
			continue
		var = [v for v in ensure_list(var) if v != unity_var]
		var_name = ''
		for i, v in enumerate(set(var)):
			var_name += v
			if var.count(v) > 1:
				var_name += f'^{var.count(v)}'
			if i == len(set(var)) - 1:
				break
			var_name += mult_sign
		if latex_mode and var_name:
			var_name = r'$\bf{' + var_name + '}$'

		if j == 0:
			formula += f'{coeff}{mult_sign if var_name else ""}{var_name} '
		else:
			if coeff > 0:
				sign = '+ '
			else:
				sign = '- '
			formula += f'{sign}{abs(coeff)}{mult_sign if var_name else ""}{var_name} '

	return formula.strip()


def molecule(molstring, mode='latex'):
	if mode == 'latex':
		ret = molstring
		for num in '0123456789':
			ret = ret.replace(num, f'_{num}')
		for sign in '+-':
			ret = ret.replace(sign, f'^{sign}')
		return ret

	if mode == 'html':
		ret = molstring
		# to take care of plus-signs used to denote reactions we have to first split 
		# the molstring into its parts.
		for part in ret.split():
			# if part is only a plus-sign we skip this part. This is only true when the plus-sign
			# is used to denote a reaction
			if part == '+':
				continue

			# parse the part
			partret = part
			# numbers should be subscript
			for num in '0123456789':
				partret = partret.replace(num, f'<sub>{num}</sub>')

			# signs should be superscript
			for sign in '+-':
				# negative charges should be denoted by em dash and not a normal dash
				partret = partret.replace(sign, f'<sup>{sign.replace("-", "—")}</sup>')

			# replace the part in the original string
			ret = ret.replace(part, partret)
		return ret


if __name__ == '__main__':
	form = polynomial([-1, 2, 3.1], ['', 'a', ('b', 'c', 'c')], latex_mode=True)
	print(form)

	print(molecule('F- + CH3Cl', 'html'))
