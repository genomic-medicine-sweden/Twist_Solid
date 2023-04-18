from jinja2 import Template
import pathlib
import time


def get_sample_name(filename):
    return pathlib.Path(filename).name.split(".")[0]


def create_report(template_filename, json_filename, show_table, tc, tc_method):
    with open(template_filename) as f:
        template = Template(source=f.read())

    with open(json_filename) as f:
        json_string = f.read()

    return template.render(
        dict(
            json=json_string,
            metadata=dict(
                date=time.strftime("%Y-%m-%d %H:%M", time.localtime()),
                sample=get_sample_name(json_filename),
                show_table=show_table,
                tc=tc,
                tc_method=tc_method,
            ),
        )
    )


def main():
    json_filename = snakemake.input.json
    template_filename = snakemake.input.template
    html_filename = snakemake.output.html

    report = create_report(
        template_filename,
        json_filename,
        snakemake.params.show_table,
        snakemake.params.tc,
        snakemake.params.tc_method,
    )

    with open(html_filename, "w") as f:
        f.write(report)


if __name__ == "__main__":
    main()
