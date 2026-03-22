class CaptionExtractor:
    def __init__(self):
        pass

    def extract(self, caption_comment):
        has_caption = caption_comment.find('caption:') != -1
        if has_caption:
            caption = caption_comment.split('caption: "')[1].split('"')[0]
        
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
        
        if has_caption:
            return f'<figcaption class="wp-element-caption">{caption}</figcaption>\n'
        else:
            return ''
