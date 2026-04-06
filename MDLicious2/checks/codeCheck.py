from MDLicious2.checks.baseCheck import BaseCheck

class CodeCheck(BaseCheck):
    def __init__(self):
        pass

    def run(self, line, warnings, errors):
        num_backticks = line.count('`')

        is_single_word = len(line.split(' ')) == 1

        if is_single_word:
            if num_backticks % 3 != 0:
                highlighted_line = self._highlight('`', line)
                errors.append(self._create_error('Code block must contain three backticks', highlighted_line))
        elif num_backticks % 6 != 0:
            highlighted_line = self._highlight('`', line)
            warnings.append(self._create_error('possible imblanaced number of backticks', highlighted_line))