from MDLicious2 import CommandLineArguments, FileProcessor, Mark2HTML, Preprocessor, CaptionMatcher
from MDLicious2 import ComponentManager, YouTube, EquationEnvironment, Code, Figure, Table
from MDLicious2 import BaseCheck, CheckManager, EquationCheck, RefCheck, CodeCheck


def main():
    # process input file
    args = CommandLineArguments()
    filereader = FileProcessor(args)

    # Build up map to to replace tags later in HTML
    caption_matcher = CaptionMatcher()
    content = caption_matcher.add_default_tag_to_equations(filereader.content)
    caption_matcher.setup_ref_map(content)

    # register customised markdown to html converters
    component_manager = ComponentManager()
    component_manager.register(YouTube(content))
    component_manager.register(EquationEnvironment(content))
    component_manager.register(Code(content))
    component_manager.register(Figure(content))
    component_manager.register(Table(content))

    # preprocess markdown file with customised converters
    # they will already convert markdown to html
    preprocessor = Preprocessor(filereader.content, component_manager.components)
    preprocessor.preprocess_content()

    # replace tags by numbers    
    content = caption_matcher.substitute(preprocessor.processed_content)

    # convert remaining (standard) markdown to html
    converter = Mark2HTML(content)
    html_content = converter.convert()

    # perform checks on input and output file
    checks = CheckManager(filereader.content, html_content, args.output)
    checks.register_markdown_checks(EquationCheck())
    checks.register_markdown_checks(CodeCheck())
    checks.register_html_checks(RefCheck())
    checks.run_checks()

    # write output file
    filereader.output(html_content)


if __name__ == '__main__':
    main()