from MDLicious2.checks.baseCheck import BaseCheck

class RefCheck(BaseCheck):
    def __init__(self):
        pass

    def run(self, line, warnings, errors):
        if '\\ref' in line:
            highlighted_line = self._highlight('\\ref', line)
            errors.append(self._create_error('ref tag found', highlighted_line))