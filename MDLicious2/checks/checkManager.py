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
        
        stderr = {}
        stderr['warnings'] = self.warnings
        stderr['errors'] = self.errors

        stderr_json = json.dumps(stderr)
        
        with open(join(self.output_directory, 'stderr.json'), 'w') as f:
            f.write(stderr_json)
