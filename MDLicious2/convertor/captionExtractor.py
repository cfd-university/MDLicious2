from enum import Enum, auto

class ComponentType(Enum):
    FIGURE = auto()
    TABLE = auto()
    CODE = auto()
    EQUATION = auto()

class CaptionExtractor:
    def __init__(self):
        self.counter = {
            ComponentType.FIGURE: 1,
            ComponentType.TABLE: 1,
            ComponentType.CODE: 1
        }

    def extract(self, caption_comment, component_type):
        has_caption = caption_comment.find('caption:') != -1
        has_tag = caption_comment.find('\\tag') != -1

        # code has a special behaviour. Check if caption and/or tag is provided and handle accordingly
        if component_type == ComponentType.CODE:
            if not has_caption and not has_tag:
                return ''

        caption = ''
        if has_caption:
            caption = caption_comment.split('caption: "')[1].split('"')[0]

        # check if caption contains a link
        if has_caption:
            caption_has_link = True
            while caption_has_link:
                caption_has_link = (caption.find('[') != -1) and (caption.find('](') != -1)
                if caption_has_link:
                    link_text = caption.split('[')[1].split('](')[0]
                    link_source = caption.split('](')[1].split(')')[0]
                    html_link = f'<a href="{link_source}" target="_blank" rel="noopener noreferrer">{link_text}</a>'
                    caption = caption.replace(f'[{link_text}]({link_source})', html_link)
                else:
                    break
        
        caption_counter = '<strong>'
        if component_type == ComponentType.FIGURE:
            figure_counter = self.counter[ComponentType.FIGURE]
            caption_counter += f'Figure {figure_counter}'
        elif component_type == ComponentType.TABLE:
            table_counter = self.counter[ComponentType.TABLE]
            caption_counter += f'Table {table_counter}'
        elif component_type == ComponentType.CODE:
            code_counter = self.counter[ComponentType.CODE]
            caption_counter += f'Listing {code_counter}'
        caption_counter += '</strong>'
        
        self.counter[component_type] += 1

        # handle all other cases
        if len(caption.strip()) > 0:
            return f'<figcaption class="wp-element-caption">{caption_counter}: {caption}</figcaption>\n'
        else:
            return f'<figcaption class="wp-element-caption">{caption_counter}</figcaption>\n'
