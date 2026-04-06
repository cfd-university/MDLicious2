from MDLicious2.checks.baseCheck import BaseCheck

class RefCheck(BaseCheck):
    def __init__(self):
        pass

    def run(self, line, warnings, errors):
        if r'\ref' in line:
            highlighted_line = self._highlight(r'\ref', line)
            errors.append(self._create_error('ref tag found', highlighted_line))

        if r'\tag' in line:
            has_correct_tag = line.find(r'\tag{eq:') != -1 or \
                              line.find(r'\tag{fig:') != -1 or \
                              line.find(r'\tag{code:') != -1 or \
                              line.find(r'\tag{tab:') != -1
            
            if not has_correct_tag:
                highlighted_line = self._highlight(r'\tag', line)
                errors.append(self._create_error(r'wrong tag found, nees to start with \ref{eq:, \ref{fig:, or \ref{tab:', highlighted_line))