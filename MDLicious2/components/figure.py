from MDLicious2.components.base import Component

class Figure(Component):
    def __init__(self, content):
        super().__init__(content)

    def match(self, index):
        line = self.content[index]
        is_figure = line.find('<!-- figure:') != -1
        has_source = False

        if is_figure:
            next_line = self.content[index + 1]
            has_source = next_line.find('![') != -1

        return is_figure and has_source

    def convert(self, index):
        # since figures take two lines, we need to increment the index by two to avoid converting the figure twice
        self.increment = 2

        # get the next two lines for processing
        line = self.content[index]
        next_line = self.content[index + 1]

        # extract meta information from image
        width = line.split('width: ')[1].split('px')[0]
        alt_text = next_line.split('![')[1].split('](')[0]
        image_source = next_line.split('](')[1].split(')')[0]

        # check if captions need to be added
        has_caption = line.find('caption:') != -1
        if has_caption:
            caption = line.split('caption: "')[1].split('"')[0]

            # check if caption contains a link
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

        figure_html = '<figure class="wp-block-image aligncenter size-large is-resized">'
        figure_html += f'<img src="{image_source}" alt="{alt_text}" class="wp-image-5550" style="width:{width}px"/>'
        if has_caption:
            figure_html += f'<figcaption class="wp-element-caption">{caption}</figcaption>'
        figure_html += '</figure>'

        return figure_html
