import json
from jinja2 import Template

def create_report(template_filename, json_filename):
    with open(template_filename) as f:
        template = Template(source=f.read())

    with open(json_filename) as f:
        json_string = f.read()

    return template.render(dict(
        json=json_string
    ))

def main():
    json_filename = snakemake.input.json
    template_filename = snakemake.input.template
    html_filename = snakemake.output.html

    report = create_report(template_filename, json_filename)

    with open(html_filename, "w") as f:
        f.write(report)

if __name__ == "__main__":
    main()
