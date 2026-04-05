from MDLicious2.checks.baseCheck import BaseCheck

class CodeCheck(BaseCheck):
    def __init__(self):
        pass

    def run(self, line, warnings, errors):
        num_backticks = line.count('`')

        if num_backticks % 6 != 0:
            highlighted_line = self._highlight('`', line)
            warnings.append(self._create_error('possible imblanaced number of backticks', highlighted_line))