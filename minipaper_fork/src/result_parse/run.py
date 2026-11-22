# main.py
from convert import markdown_to_html
from pathlib import Path


def main():
    parsed_dir = Path("/home/kido/avinasi/clockbase/data/parsed_output")
    for gse in parsed_dir.glob("*"):
        gse_id = gse.name
        md_path = gse / "parsed_report.md"
        html_path = gse / "parsed_report.html"
        markdown_content = md_path.read_text(encoding="utf-8")
        try:
            html_output = markdown_to_html(markdown_content, gse_id)
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_output)
        except Exception as e:
            raise e
            print(f"Error processing {md_path}: {e}")
            with open('error_html.log', 'a') as f:
                f.write(f"{md_path}: {e}\n")


if __name__ == "__main__":
    main()