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
		for part in ret.split():
			if part in '+-':
				continue

			partret = part
			for num in '0123456789':
				partret = partret.replace(num, f'<sub>{num}</sub>')
			for sign in '+-':
				# if sign == '-'
				partret = partret.replace(sign, f'<sup>{sign.replace("-", "—")}</sup>')
			ret = ret.replace(part, partret)
		return ret


if __name__ == '__main__':
	form = polynomial([-1, 2, 3.1], ['', 'a', ('b', 'c', 'c')], latex_mode=True)
	print(form)

	print(molecule('F- + CH3Cl', 'html'))