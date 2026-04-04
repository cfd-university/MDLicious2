from MDLicious2 import CommandLineArguments, FileProcessor, Mark2HTML, Preprocessor, CaptionMatcher
from MDLicious2 import ComponentManager, YouTube, EquationEnvironment, Code, Figure, Table
from MDLicious2 import BaseCheck, CheckManager, EquationCheck


def main():
    # process input file
    args = CommandLineArguments()
    filereader = FileProcessor(args)

    # perform checks on input file
    checks = CheckManager(filereader.content, args.output)
    checks.register(EquationCheck())
    checks.run_checks()

    # Replace \ref{} call with actual captions
    caption_matcher = CaptionMatcher(filereader.content)
    content = caption_matcher.substitute()

    # register customised markdown to html converters
    component_manager = ComponentManager()
    component_manager.register(YouTube(content))
    component_manager.register(EquationEnvironment(content))
    component_manager.register(Code(content))
    component_manager.register(Figure(content))
    component_manager.register(Table(content))

    # preprocess markdown file with customised converters
    # they will already convert markdown to html
    preprocessor = Preprocessor(content, component_manager.components)
    preprocessor.preprocess_content()
        
    # convert remaining (standard) markdown to html
    converter = Mark2HTML(preprocessor.processed_content)
    html_content = converter.convert()    

    # write output file
    filereader.output(html_content)


if __name__ == '__main__':
    main()