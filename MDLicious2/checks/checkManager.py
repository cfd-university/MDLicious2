from os.path import join
import json

class CheckManager:
    def __init__(self, content, output_directory):
        self.content = content
        self.output_directory = output_directory

        self.checks = []
        self.warnings = []
        self.errors = []

    def register(self, check):
        self.checks.append(check)

    def run_checks(self):
        for check in self.checks:
            for line in self.content:
                check.run(line, self.warnings, self.errors)
        
        warning_json = json.dumps(self.warnings)
        error_json = json.dumps(self.errors)
        
        with open(join(self.output_directory, 'warnings.json'), 'w') as f:
            f.write(warning_json)
        
        with open(join(self.output_directory, 'errors.json'), 'w') as f:
            f.write(error_json)