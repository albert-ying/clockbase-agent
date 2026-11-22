import markdown_it
from mdit_py_plugins.front_matter import front_matter_plugin
from mdit_py_plugins.dollarmath import dollarmath_plugin

def markdown_to_html(markdown_text: str, gse_id) -> str:
    """
    Convert academic-style Markdown text into a self-contained, well-formatted HTML file.

    Args:
    markdown_text (str): A string containing Markdown text.

    Returns:
    str: A complete HTML-formatted string.
    """
    
    # Initialize markdown-it and enable plugins
    # - dollarmath_plugin: for parsing $...$ and $$...$$ LaTeX math expressions
    # - front_matter_plugin: for parsing YAML front matter (e.g. --- title: ---)
    md = (
        markdown_it.MarkdownIt("default", {"html": True, "linkify": True, "typographer": True})
        .use(front_matter_plugin)
        .use(dollarmath_plugin)
    )

    # Render Markdown to HTML
    html_body = md.render(markdown_text)

    # Extract title from front matter (if present)
    # Default title
    title = gse_id
    
    # Try to get front matter from the environment
    env = {}
    md.parse(markdown_text, env)
    if "frontMatter" in env:
        title = env["frontMatter"].get("title", gse_id)

    # HTML template including CSS and JavaScript
    # CSS: for overall aesthetics, typography, and element styling
    # KaTeX: for rendering LaTeX math expressions
    # highlight.js: for code syntax highlighting
    html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css" integrity="sha384-n8MVd4RsNIU0KOVEMeaOVTfZhFtYgiEjZFSln3eLjGqKCMOomLQB0Eg9joPlsOFQ" crossorigin="anonymous">
    
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/highlight.js@11.9.0/styles/atom-one-light.min.css">

    <style>
        /* Google Fonts for elegant typography */
        @import url('https://fonts.googleapis.com/css2?family=Source+Sans+Pro:wght@400;600;700&family=Source+Serif+Pro:ital,wght@0,400;0,600;1,400&display=swap');

        body {{
            font-family: 'Source Serif Pro', serif;
            line-height: 1.7;
            color: #333;
            background-color: #fdfdfd;
            margin: 0;
            padding: 20px;
            display: flex;
            justify-content: center;
        }}

        main {{
            max-width: 800px;
            width: 100%;
            padding: 2em;
            background-color: #fff;
            box-shadow: 0 0 20px rgba(0,0,0,0.05);
            border-radius: 8px;
        }}

        h1, h2, h3, h4, h5, h6 {{
            font-family: 'Source Sans Pro', sans-serif;
            color: #1a1a1a;
            font-weight: 600;
            margin-top: 1.5em;
            margin-bottom: 0.5em;
            line-height: 1.2;
        }}

        h1 {{
            font-size: 2.5em;
            border-bottom: 1px solid #eee;
            padding-bottom: 0.3em;
            margin-top: 0;
            line-height: 1.1;
        }}

        h2 {{
            font-size: 2em;
            border-bottom: 1px solid #eee;
            padding-bottom: 0.3em;
        }}

        p {{
            margin: 1em 0;
        }}

        a {{
            color: #007bff;
            text-decoration: none;
            transition: color 0.2s;
        }}

        a:hover {{
            color: #0056b3;
            text-decoration: underline;
        }}

        img {{
            max-width: 100%;
            height: auto;
            display: block;
            margin: 1.5em auto;
            border-radius: 4px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }}

        blockquote {{
            margin: 1.5em 0;
            padding: 0.5em 1.5em;
            border-left: 5px solid #007bff;
            background-color: #f8f9fa;
            color: #555;
        }}

        code {{
            font-family: 'Menlo', 'Consolas', monospace;
            background-color: #f1f1f1;
            padding: 0.2em 0.4em;
            border-radius: 3px;
            font-size: 0.9em;
        }}

        pre {{
            background-color: #f8f9fa;
            border: 1px solid #e1e1e1;
            padding: 1em;
            border-radius: 4px;
            overflow-x: auto;
        }}

        pre code {{
            background-color: transparent;
            padding: 0;
            font-size: 0.9em;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1.5em 0;
            box-shadow: 0 0 10px rgba(0,0,0,0.05);
        }}

        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}

        th {{
            background-color: #f8f9fa;
            font-family: 'Source Sans Pro', sans-serif;
            font-weight: 600;
        }}

        /* Katex display math centering and scroll */
        .katex-display {{
            display: block;
            text-align: center;
            margin: 1em 0;
            overflow-x: auto;
            overflow-y: hidden;
            padding: 0.5em 0;
        }}

    </style>
</head>
<body>
    <main>
        {html_body}
    </main>

    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js" integrity="sha384-XjKyOOlGwcjNTAIQHIpgOno0Hl1YQqzUOEleOLALmuqehneUG+vnGctmUb0ZY0l8" crossorigin="anonymous"></script>
    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/contrib/auto-render.min.js" integrity="sha384-+VBxd3r6XgURycqtZ117nYw44OOcIax56Z4dCRWbxyPt0Koah1uHoK0o4+/RRE05" crossorigin="anonymous"></script>
    
    <script src="https://cdn.jsdelivr.net/npm/highlight.js@11.9.0/lib/common.min.js"></script>

    <script>
        // Initialize KaTeX auto-rendering
        document.addEventListener("DOMContentLoaded", function() {{
            renderMathInElement(document.body, {{
                delimiters: [
                    {{left: '$$', right: '$$', display: true}},
                    {{left: '$', right: '$', display: false}},
                    {{left: '\\\\(', right: '\\\\)', display: false}},
                    {{left: '\\\\[', right: '\\\\]', display: true}}
                ],
                throwOnError : false
            }});
            
            // Initialize highlight.js
            hljs.highlightAll();
        }});
    </script>
</body>
</html>
"""
    return html_template