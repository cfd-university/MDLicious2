from os.path import join
import json

class CheckManager:
    def __init__(self, markdown, html, output_directory):
        self.markdown = markdown
        self.html = html
        self.output_directory = output_directory

        self.pre_checks = []
        self.post_checks = []
        self.warnings = []
        self.errors = []

    def register_markdown_checks(self, check):
        self.pre_checks.append(check)

    def register_html_checks(self, check):
        self.post_checks.append(check)

    def run_checks(self):
        # check that markdown input does not contain errors
        for check in self.pre_checks:
            for line in self.markdown:
                check.run(line, self.warnings, self.errors)
        
        # check that html output does not contain errors
        for check in self.post_checks:
            for line in self.markdown:
                check.run(line, self.warnings, self.errors)

        stderr = {}
        stderr['warnings'] = self.warnings
        stderr['errors'] = self.errors

        stderr_json = json.dumps(stderr, indent=2)
        
        with open(join(self.output_directory, 'stderr.json'), 'w') as f:
            f.write(stderr_json)
