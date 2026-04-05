from MDLicious2.checks.baseCheck import BaseCheck

class EquationCheck(BaseCheck):
    def __init__(self):
        pass

    def run(self, line, warnings, errors):
        num_dollar_signs = line.count('$')
        is_even = num_dollar_signs % 2 == 0

        if not is_even:
            highlighted_line = self._highlight('$', line)
            warnings.append(self._create_warning('uneven number of dollar signs', highlighted_line))

        if len(line) > 2:
            double_dollar_signs = line.count('$$')

            if double_dollar_signs > 0:
                highlighted_line = self._highlight('$$', line)
                errors.append(self._create_error('double dollar signs', highlighted_line))