from MDLicious2.checks.baseCheck import BaseCheck

class EquationCheck(BaseCheck):
    def __init__(self):
        pass

    def run(self, line, warnings, errors):
        num_dollar_signs = line.count('$')
        is_even = num_dollar_signs % 2 == 0

        if not is_even:
            highlighted_line = line.replace('$', '<strong>$</strong>')
            warnings.append(f'uneven number of dollar signs in line: {highlighted_line}')

        if len(line) > 2:
            double_dollar_signs = line.count('$$')

            if double_dollar_signs > 0:
                highlighted_line = line.replace('$$', '<strong>$$</strong>')
                errors.append(f'double dollar signs in line: {highlighted_line}')